"""Define the base data containing classes for precisionguide system

These are designed to be portable across any platform, since they do
not depend on the specifics of the telescsope control system.  Thus,
data can be recorded using this code on, e.g., a Windows system and
uploaded to any other platform for detailed analysis using this same
code.

"""
import numpy as np
from scipy.signal import convolve2d
from scipy.ndimage.measurements import center_of_mass

from astropy.io import fits
from astropy.nddata import CCDData, fits_ccddata_reader
from astropy import units as u
from astropy.time import Time
from astropy.convolution import Gaussian2DKernel

from IoIO.ccdmultipipe import fallback_unit_ccddata_reader

_NotFound = object()

class pgproperty(property):
    """Caching property decorator with auto setting and None reset

    . myprop = None resets system and forces the getter to run

    . when getter runs and returns a non-None value, that value is
    cached and returned instead of running the getter

    . myprop = 'some_value' overrides the getter -- 'some_value' is
    cached and returned on subsequent queries of myprop

    . No setter is required to set the property, but if one is
    provided, its return value is used as the cached value that
    overrides the getter (i.e., not knowledge of the internal workings
    of the system is requried to make it work)

    . Deleters can be specified to do something before the cache
    reference is deleted.  Deletion is not permanent -- the next time
    myprop is queried, the getter is run

    . shape_check class variable: None, 0, or tuple

        This class variable affects the treatment of the value
        returned by the setter or if no setter is
        present, value provided by user.  
        . None: the value is cached without modification
        . 0: value is converted to a numpy array with np.asarray()
        . tuple: shape of resulting np.array must match tuple or a
        ValueError is raised

    Inspired by `astropy.utils.lazyproperty` and `property_manager
    <https://github.com/xolox/python-property-manager>`

    .NOTE: This cache operates on the top-level property.  In a
    complicated multi-inheritance system with inter-connected property
    like `pgdata`, some time savings in calculation of quantities
    could be achieved by caching property at each inheritance level
    individually, say by making the `key` a string that includes the
    property name, `self.__class__.__name__` and also the module name
    just to make sure it is unique.  Then a None reset would only
    reset the levels from the top down to the MRO at which the
    property was set to None.  This guarantees the children get a
    freshly calculated value of the affected properties and
    acknowledges that the parents don't care that the change was made,
    since they never depended on the relationship between the property
    in the first place.  Or if they did, they should be super()ed into
    the calculation and the most senior setter concerned about the
    interrelationship would also be calling for the reset of the
    affected other property.  The side-affect of caching all of the
    intermediate results would be more memory use, since all of the
    intermediate property would have long-lived references.  For the
    pgdata coordiate-oriented stuff, this is not a big deal, but it is
    not generally appropriate

    """
    
    shape_check = None
    def __init__(self, fget,
                 fset=None,
                 fdel=None,
                 doc=None):
        super().__init__(fget, fset, fdel, doc)
        self._key = self.fget.__name__
        self._none_proxy = 'None proxy'

    def npval(self, val):
        """Turn val into np.array if shape_check is non-None.  Checks
        shape of resulting array if shape_check is a tuple"""
        if self.shape_check is None:
            return val
        val = np.asarray(val)
        if (self.shape_check != 0
            and val.shape != self.shape_check):
            raise ValueError(f'value "{val}" is not a sensible {self.shape_check}-dimensional coordinate')
        return val

    def __get__(self, obj, owner=None):
        try:
            obj_dict = obj.__dict__
            val = obj_dict.get(self._key, _NotFound)
            if val is self._none_proxy:
                return None
            if val is _NotFound or val is None:
                val = self.fget(obj)
                if val is None:
                    obj_dict[self._key] = self._none_proxy
                    return None
                val = self.npval(val)
                obj_dict[self._key] = val
            return val
        except AttributeError:
            if obj is None:
                return self
            raise

    def __set__(self, obj, val):
        obj_dict = obj.__dict__
        if self.fset:
            val = self.fset(obj, val)
        if val is not None:
            val = self.npval(val)
        obj_dict[self._key] = val

    def __delete__(self, obj):
        if self.fdel:
            self.fdel(obj)
        obj.__dict__.pop(self._key, None)    # Delete if present


class pgcoordproperty(pgproperty):
    shape_check=(2,)

class PGCenter():
    """Base class for containing object center and desired center


    Parameters
    ----------
    obj_center : tuple or array
        Center of object, (y, x) in *unbinned* coordinates referenced to the
        *origin* of the CCD.  Note Python C array indexing into 2D array

    desired_center : tuple-like
        Desired center of object, (y, x) in *unbinned* coordinates
        referenced to the *origin* of the CCD.  Note Python C array
        indexing into 2D array

    quality : int
        0 -- 10 ranking of quality of quality of ``obj_center``
        determination, with 0 = very bad, 10 = very good

    tmid : `~astropy.time.Time`
        Time at midpoint of observation

    """
    def __init__(self,
                 obj_center=None,
                 desired_center=None,
                 quality=None,
                 tmid=None):
        self.obj_center = obj_center
        self.desired_center = desired_center
        self.quality = quality
        self.tmid = tmid

    # Our decorator does everything we need :-)
    @pgcoordproperty
    def obj_center(self):
        pass
        
    @pgcoordproperty
    def desired_center(self):
        pass

    @pgproperty
    def quality(self):
        pass
        
    @quality.setter
    def quality(self, value):
        """Unset quality translates to 0 quality"""
        if value is None:
            value = 0
        if not isinstance(value, int) or value < 0 or value > 10:
            raise ValueError('quality must be an integer value from 0 to 10')
        else:
            return value
        
    @pgcoordproperty
    def tmid(self):
        pass

class FBUCCDData(CCDData):
    """
    Add ``fallback_unit`` capability to :meth:`CCDData.read() <~astropy.nddata.CCDData.read>`

    Paramters
    ---------
    filename : str
        Name of FITS file

    unit : `~astropy.units.Unit`, optional
        Units of the image data.   If ``None,`` the FITS header
        ``BUNIT`` keyword will be queried to set the unit.  If that
        is not found or an error is raised, `fallback_unit` will be
        used.
        Default is ``None``.

    fallback_unit : `~astropy.units.Unit`, optional
        Units to be used for the image data if `unit` is not provided
        and no valid ``BUNIT`` value is found in the FITS header.
        Default is ``None``.

    kwargs :
        Keywords to pass to :meth:`CCDData.read() <~astropy.nddata.CCDData.read>`

    """
    def __init__(self, *args,
                 filename=None,
                 fallback_unit=None,
                 **kwargs):
        if filename is not None:
            ccd = fallback_unit_ccddata_reader(filename, *args, 
                                               fallback_unit=fallback_unit,
                                               **kwargs)
            self.__dict__.update(ccd.__dict__)
        else:
            super().__init__(*args, **kwargs)

    @classmethod
    def read(cls, filename, *args,
             fallback_unit=None,
             **kwargs):
        return cls(*args,
                   filename=filename, 
                   fallback_unit=fallback_unit,
                   **kwargs)
        #ccd = fallback_unit_ccddata_reader(filename, *args, 
        #                                   fallback_unit=fallback_unit,
        #                                   **kwargs)
        #return ccd

class PGData(CCDData):
    """Base class for image data in the `precisionguide` system

    This class stores CCD data using `~astropy.nddata.CCDData` as a
    superclass and adds properties and methods that calculate/store
    four quantities: :prop:`obj_center`, :prop:`desired_center`,
    :prop:`quality`, and :prop:`tmid`.  These four quantities are
    intended to be returned in a :class:`PGCenter` object for
    subsequent lightweight storage and use in the precisionguide
    system.  Because precisionguide controls the absolute position of
    an object (or FOV center) on the CCD, :prop:`desired_center` and
    :prop:`obj_center` always read in *unbinned* pixel values
    referenced to the origin of the CCD itself.  Thus, the image input
    to :class:`PGData` must include both the image FITS header and
    image array(s)

    Parameters
    ----------

    """

    def __init__(self,
                 *args,
                 unit=None,
                 raw_unit='adu',
                 obj_center=None,
                 desired_center=None,
                 quality=0,
                 date_obs_key='DATE-OBS',
                 exptime_key='EXPTIME',
                 darktime_key='DARKTIME',
                 **kwargs):
        if unit is None:
            unit = raw_unit
        super().__init__(*args, unit=unit, **kwargs)
        self.obj_center = obj_center
        self.desired_center = desired_center
        self.quality = quality
        self.date_obs_key = date_obs_key
        self.exptime_key = exptime_key
        self.darktime_key = darktime_key

    @classmethod
    def read(cls, filename, *args, 
             raw_unit='adu',
             obj_center=None,
             desired_center=None,
             quality=0,
             date_obs_key='DATE-OBS',
             exptime_key='EXPTIME',
             darktime_key='DARKTIME',
             **kwargs):
        #ccd = fallback_unit_ccddata_reader(filename, *args,
        ccd = FBUCCDData.read(filename, *args, 
                              fallback_unit=raw_unit,
                              **kwargs)
        # Make a vestigial pgd
        pgd = cls(ccd.data,
                  unit=ccd.unit,
                  obj_center=obj_center,
                  desired_center=desired_center,
                  quality=quality,
                  date_obs_key=date_obs_key,
                  exptime_key=exptime_key,
                  darktime_key=darktime_key)
        # Merge in ccd property
        pgd.__dict__.update(ccd.__dict__)
        return pgd

    @pgcoordproperty
    def obj_center(self):
        # Here is our fancy calculation that takes a long time.
        # DESIGN NOTE: obj_center and desired_center use 
        # the pgproperty system so that they can be set to an
        # arbitrary value at any time.  Recalculation is triggered by 
        # and stored as an np.array without recalculation.  In
        # other words, we calculate just once but we enable 
        obj_center = np.asarray((-99, -99))
        self.quality = 0
        return obj_center

    @pgcoordproperty
    def desired_center(self):
        # Here is where the complicated desired_center calculation is
        # done:
        return np.asarray(self.data.shape)/2

    @pgproperty
    def quality(self):
        self.obj_center

    @pgproperty
    def tmid(self):
        tmid_str = self.meta.get('tmid')
        if tmid_str is not None:
            return Time(tmid_str, format='fits')
        try:
            exptime = self.meta.get(self.darktime_key.lower()) 
            if exptime is None:
                exptime = self.meta[self.exptime_key.lower()]
            exptime *= u.s
            dateobs_str = self.meta[self.date_obs_key.lower()] 
            return Time(dateobs_str, format='fits') + exptime/2
        except:
            log.error(f'Cannot read sufficient combination of '
                      '{self.date_obs_key}, {self.darktime_key} '
                      'and/or {self.exptime_key} keywords from FITS '
                      'header to establish tmid')
            return None

    @tmid.setter
    def tmid(self, val):
        if not isinstance(val, Time):
            raise ValueError('tmid must be of type `~astropy.time.Time`')

    @property
    def binning(self):
        """Image binning in Y,X order"""
        # Note, this needs to be overwritten with the actual FITS
        # binning keywords used (which are unfortunately not in any
        # FITS standard).  E.g.:
        #binning = np.asarray((self.meta['YBINNING'],
        #                      self.meta['XBINNING']))
        binning = (1,1)
        return np.asarray(binning)
        
    @property
    def subframe_origin(self):
        """Subframe origin in *unbinned* pixels with full CCD origin = (0,0).  Y,X order"""
        #subframe_origin = np.asarray((self.meta['YORGSUBF'],
        #                              self.meta['XORGSUBF']))
        #subframe_origin *= self.binning
        subframe_origin = (0,0)
        return np.asarray(subframe_origin)

    def unbinned(self, coords):
        """Returns coords referenced to full CCD given internally stored binning/subim info"""
        coords = np.asarray(coords)
        return np.asarray(self.binning * coords + self.subframe_origin)

    def binned(self, coords):
        """Assuming coords are referenced to full CCD, return location in binned coordinates relative to the subframe origin"""
        coords = np.asarray(coords)
        return np.asarray((coords - self.subframe_origin) / self.binning)
        
    def im_unbinned(self, a):
        """Returns an unbinned version of a.  a must be same shape as data
        """
        assert a.shape == self.data.shape
        # Don't bother if we are already unbinned
        if np.sum(self.binning) == 2:
            return a
        newshape = self.binning * a.shape
        # From http://scipy-cookbook.readthedocs.io/items/Rebinning.html
        assert len(a.shape) == len(newshape)

        slices = [ slice(0,old, float(old)/new)
                   for old,new in zip(a.shape,newshape) ]
        coordinates = np.mgrid[slices]
        indices = coordinates.astype('i')   #choose the biggest smaller integer index
        unbinned = a[tuple(indices)]
        # Check to see if we need to make a larger array into which to
        # plop unbinned array
        if np.sum(self.subframe_origin) > 0:
            # Note subframe origin reads in binned pixels
            origin = self.unbinned(self.subframe_origin)
            full_unbinned = np.zeros(origin + np.asarray(unbinned.shape))
            full_unbinned[origin[0]:, origin[1]:] = unbinned
            unbinned = full_unbinned
        return unbinned

    @pgproperty
    def data_unbinned(self):
        """Returns an unbinned version of data
        """
        return self.im_unbinned(self.data)

    @property
    def pgcenter(self):
        return PGCenter(self.unbinned(self.obj_center),
                        self.unbinned(self.desired_center),
                        self.quality,
                        self.tmid)

    def _card_write(self):
        # Note pythonic y, x coordinate ordering
        self.meta['DES_CR0'] = (self.desired_center[1], 'Desired center X')
        self.meta['DES_CR1'] = (self.desired_center[0], 'Desired center Y')
        self.meta['OBJ_CR0'] = (self.obj_center[1], 'Object center X')
        self.meta['OBJ_CR1'] = (self.obj_center[0], 'Object center Y')
        self.meta['QUALITY'] = (self.quality, 'Quality on 0-10 scale of center determination')
        self.meta['QUALITY'] = (self.quality, 'Quality on 0-10 scale of center determination')
        self.meta.insert('DATE-OBS',
                         ('TMID', self.tmid.value,
                          'midpoint of exposure, UT'),
                         after=True)

    def write(self, *args, **kwargs):
        self._card_write()
        super().write(*args, **kwargs)

class NoCenterPGD(PGData):
    """Sets :prop:`obj_center` to invalid value and :prop:`quality` to 0`"""
    pass

class CenteredPGD(PGData):
    """Sets :prop:`obj_center` to :prop:`desired_center`"""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._old_desired_center = None
    

    @pgcoordproperty
    def obj_center(self):
        self.quality = 10
        return self.desired_center

    @pgcoordproperty
    def desired_center(self):
        return super().desired_center

    @desired_center.setter
    def desired_center(self, value):
        # Reset obj_center for proper recalculation
        del self.obj_center
        return value

class OffsetPGD(PGData):
    """Offsets :prop:`obj_center` by :prop:`center_offset`  Note:
    proper object that determins center needs to be inherited before
    this for sensible center: e.g.: 
    class MyPGD(PGDOffset, PGDCentered): pass"""

    def __init__(self, *args,
                 center_offset=None,
                 **kwargs):
        super().__init__(*args, **kwargs)
        # Center offset is used in the calculation of a
        # pgcoordproperty such that we want to delete that
        # property for recalculation in the event center_offset
        # changes.  Particularly considering we end up in np.arrays
        # and the += operator increments the object numpy object
        # itself (e.g. pgd.center_offset += (10,10)), this is just a
        # little too complex for the pgcoordproperty system to
        # handle, so do it the olde fashionede waye

        if center_offset is None:
            center_offset = (0,0)
        self._old_center_offset = center_offset
        self.center_offset = center_offset
        #self._center_offset = np.asarray(center_offset)

    @pgcoordproperty
    def _old_center_offset(self):
        pass        

    @pgcoordproperty
    def center_offset(self):
        pass

    @center_offset.setter
    def center_offset(self, value):
        if value is None:
            value = (0,0)
        value = np.asarray(value)
        if not np.all(self._old_center_offset == value):
            del self.obj_center
            # This copy is key for avoiding subtle assign by reference bugs
            self._old_center_offset = value.copy()
            self.obj_center = None
        return value

    @pgcoordproperty
    def obj_center(self):
        return super().obj_center + self.center_offset

class MaxPGD(PGData):
    @pgcoordproperty
    def obj_center(self):
        self.quality = 6
        obj_center = np.unravel_index(np.argmax(self), self.shape)
        return np.asarray(obj_center)
    
class BrightestPGD(PGData):
    def __init__(self, *args,
                 seeing=2,
                 **kwargs):
        super().__init__(*args, **kwargs)
        self.seeing = seeing

    @pgcoordproperty
    def obj_center(self):
        # https://keflavich-astropy.readthedocs.io/en/convolve_fft_profiling/convolution/performance.html
        # suggests that astropy convolv is only needed when their are NaNs

        # Still working on this
        return super().obj_center + self.center_offset
    
class CenterOfMassPGD(PGData):
    """Use simple center-of-mass calculation to find center of
    brightest object.  Works best if pixels not likely to be part of
    the source are set to zero."""

    @pgcoordproperty
    def obj_center(self):
        # Still working on this
        return super().obj_center + self.center_offset

class MaxImPGD(PGData):
    """Handles MaxIM DL :prop:`binning` and :prop:`subframe_origin`"""

    # It is important to just read these once, since we may modify
    # them in the object before write
    @pgcoordproperty
    def binning(self):
        """Image binning in Y,X order"""
        binning = (self.meta['YBINNING'],
                   self.meta['XBINNING'])
        return binning
        
    @pgcoordproperty
    def subframe_origin(self):
        """Subframe origin in *unbinned* pixels with full CCD origin = (0,0).  Y,X order"""
        subframe_origin = (self.meta['YORGSUBF'],
                           self.meta['XORGSUBF'])
        subframe_origin = np.asarray(subframe_origin) * self.binning
        return subframe_origin

    def _card_write(self):
        # Note pythonic y, x coordinate ordering
        self.meta['XBINNING'] = self.binning[1]
        self.meta['YBINNING'] = self.binning[0]
        self.meta['XORGSUBF'] = self.subframe_origin[1]
        self.meta['YORGSUBF'] = self.subframe_origin[0]
        super()._card_write()

    def write(self, *args, **kwargs):
        self._card_write()
        super().write(*args, **kwargs)

#pgc = PGCenter()
#pgc.obj_center = (1,1)
#print(pgc.obj_center)
##pgc.obj_center = 1
##pgc = PGCenter((1,2), (3,4))
#pgc = PGCenter([1,2], (3,4))
#print(pgc.obj_center, pgc.desired_center)

#pgd = PGData()
fname = '/data/io/IoIO/raw/2020-07-15/HD87696-0016_Na_off.fit'
#
#class MyPGD(OffsetPGD, CenteredPGD, PGData):
#    pass

#pgd = PGData.read(fname)
#print(pgd.obj_center, pgd.desired_center)
#pgd.obj_center = (1,1)
#print(pgd.obj_center, pgd.desired_center)
#pgd = PGData.read(fname, obj_center=(2,2))
#print(pgd.obj_center, pgd.desired_center)
#pgd.desired_center = (1,1)
#print(pgd.obj_center, pgd.desired_center)
#del pgd.desired_center
#print(pgd.obj_center, pgd.desired_center)
#pgd.obj_center = None
#print(pgd.obj_center, pgd.desired_center)
#pgd = MyPGD.read(fname)
#print(pgd.obj_center, pgd.desired_center, pgd.center_offset)
#pgd.center_offset += np.asarray((10, -10))
#pgd.center_offset += (10, -10)
#print(pgd.obj_center, pgd.desired_center, pgd.center_offset)
#pgd.center_offset = None
#pgd.desired_center = (3,3)
#pgd.center_offset += (10, -30)
###pgd.center_offset = (10, -10)
#print(pgd.obj_center, pgd.desired_center, pgd.center_offset)
#pgd.obj_center = (0,0)
#print(pgd.obj_center, pgd.desired_center, pgd.center_offset)
#pgd.obj_center = None
##del pgd.obj_center
##pgd.calculated_center = (0,0)
#print(pgd.obj_center, pgd.desired_center, pgd.center_offset)
#pgd.center_offset = pgd.center_offset + (10, -10)
#print(pgd.obj_center, pgd.desired_center, pgd.center_offset)
#pgd.center_offset += (10, -30)
#print(pgd.obj_center, pgd.desired_center, pgd.center_offset)

#pgd.obj_center = (1,1)
#pgd.desired_center = (1,1)
#print(pgd.obj_center, pgd.desired_center)
#del pgd.desired_center
#print(pgd.obj_center, pgd.desired_center)
#pgd.center_offset = (20. -10)
#print(pgd.obj_center, pgd.desired_center, pgd.center_offset)

#pgd = MaxPGD.read(fname)
#print(pgd.obj_center)
class MyPGD(OffsetPGD, MaxPGD):
    pass
#pgd = MyPGD.read(fname, center_offset=(-10,20))
#print(pgd.obj_center, pgd.desired_center, pgd.center_offset)

bname = '/data/io/IoIO/reduced/Calibration/2020-07-07_ccdT_-10.3_bias_combined.fits'
bpgd = MyPGD.read(bname)
epgd = MyPGD.read(bname, unit='adu')
#ccd = FBUCCDData.read(bname)

pgd = MyPGD.read(fname)


#print(pgd.desired_center)

#
#ccd = CCDData.read(fname, unit='adu')
#
#rpgd = PGData(ccd.data, meta=ccd.meta)
#print(pgd.desired_center)
#
#pgd = PGData.read(fname)
#print(pgd.desired_center)
#
#pgd = PGData.read(fname, desired_center=(2,2))
#print(pgd.desired_center)
#
#dspgd = PGData.read(fname, desired_center=(100,100))
#print(dspgd.obj_center, dspgd.desired_center)
#
#dspgd.write('/tmp/test.fits', overwrite=True)
#
#pgdc = PGDCentered.read(fname)
#print(pgdc.obj_center, pgdc.desired_center)
#
#pgdo = PGDOffset.read(fname)
#print(pgdo.obj_center, pgdo.desired_center)
#
#pgdo = PGDOffset(pgd.data, meta=pgd.meta, center_offset=(20,10))
#print(pgdo.obj_center, pgdo.desired_center)
#
#pgdo = PGDOffset.read(fname, center_offset=(20,10))
#print(pgdo.obj_center, pgdo.desired_center)
#
#class MyPGD(PGDOffset, PGDCentered):
#    pass
#
#mpgd = MyPGD.read(fname, center_offset=(20,10))
#print(mpgd.obj_center, mpgd.desired_center)
#
#print('done')

bname = '/data/io/IoIO/reduced/Calibration/2020-07-07_ccdT_-10.3_bias_combined.fits'
#ccd = CCDData.read(bname)

fname1 = '/data/Mercury/raw/2020-05-27/Mercury-0005_Na-on.fit'
#ccd = CCDData.read(fname1)
ccd = FBUCCDData.read(fname1, fallback_unit='adu')
#ccd = FBUCCDData.read(fname1, fallback_unit='aduu')

ccd = FBUCCDData.read(fname1, unit='electron')
