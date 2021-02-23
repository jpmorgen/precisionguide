"""Define the base data containing classes for precisionguide system

These are designed to be portable across any platform, since they do
not depend on the specifics of the telescsope control system.  Thus,
data can be recorded using this code on, e.g., a Windows system and
uploaded to any other platform for detailed analysis using this same
code.

"""

import numpy as np

from astropy.nddata import CCDData
from astropy import units as u
from astropy.time import Time

class PGCenter():
    """Base class for containing object center and desired center


    Parameters
    ----------
    obj_center : 
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

    @property
    def obj_center(self):
        return self._obj_center
        
    @obj_center.setter
    def obj_center(self, value):
        if value is None:
            self._obj_center = None
        else:
            self._obj_center = np.asarray(value)

    @property
    def desired_center(self):
        if self._desired_center is not None:
            return self._desired_center
    
    @desired_center.setter
    def desired_center(self, value):
        if value is None:
            self._desired_center = None
        else:
            self._desired_center = np.asarray(value)

    @property
    def quality(self):
        if self._quality is not None:
            return self._quality
        
    @quality.setter
    def quality(self, value):
        if value is None:
            self._quality = None
            return self._quality
        if not isinstance(value, int) or value < 0 or value > 10:
            raise ValueError('quality must be an integer value from 0 to 10')
        else:
            self._quality = value
        
    @property
    def tmid(self):
        if self._tmid is not None:
            return self._tmid
        
    @tmid.setter
    def tmid(self, value):
        if value is None:
            self._tmid = None
        else:
            self._tmid = np.asarray(value)


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
                 from_read=False,
                 unit=None,
                 raw_unit='adu',
                 desired_center=None,
                 date_obs_key='DATE-OBS',
                 exptime_key='EXPTIME',
                 darktime_key='DARKTIME',
                 **kwargs):
        if not from_read:
            if unit is None:
                unit = raw_unit
            super().__init__(*args, unit=unit, **kwargs)
        self._data_unbinned = None
        self._desired_center = None
        self._obj_center = None
        self.desired_center = desired_center
        self.date_obs_key = date_obs_key
        self.exptime_key = exptime_key
        self.darktime_key = darktime_key

    @classmethod
    def read(cls, *args,
             unit=None,
             raw_unit='adu',
             desired_center=None,
             date_obs_key='DATE-OBS',
             exptime_key='EXPTIME',
             darktime_key='DARKTIME',
             **kwargs):
        if unit is None:
            unit = raw_unit
        pgd = super(PGData, cls).read(*args, unit=unit, **kwargs)
        pgd.__init__(from_read=True,
                     desired_center=desired_center,
                     date_obs_key=date_obs_key,
                     exptime_key=exptime_key,
                     darktime_key=darktime_key)
        return pgd

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

    @property
    def data_unbinned(self):
        """Returns an unbinned version of data
        """
        if self._data_unbinned is not None:
            return self._data_unbinned
        self._data_unbinned = self.im_unbinned(self.data)
        return self._data_unbinned

    @property
    def obj_center(self):
        if self._obj_center is not None:
            return self._obj_center
        self._obj_center = (-99, -99)
        self.quality = 0
        return self._obj_center

    @obj_center.setter
    def obj_center(self, value):
        if value is None:
            self._obj_center = None
        else:
            self._obj_center = np.asarray(value)

    @property
    def desired_center(self):
        if self._desired_center is not None:
            return self._desired_center
        # Here is where the complicated desired_center calculation is
        # done:
        self._desired_center = np.asarray(self.data.shape)/2
        return self._desired_center

    @desired_center.setter
    def desired_center(self, value=None):
        if value is None:
            self._desired_center = None
        else:
            self._desired_center = np.asarray(value)

    @property
    def tmid(self):
        tmid = self.meta.get('tmid')
        if self.meta.get('tmid') is not None:
            return tmid
        try:
            exptime = self.meta.get(self.darktime_key.lower()) 
            if exptime is None:
                exptime = self.meta[self.exptime_key.lower()]
            exptime *= u.s
            tmid = (Time(self.meta[self.date_obs_key.lower()],
                         format='fits')
                    + exptime/2)
        except:
            log.warning(f'Cannot read {self.darktime_key} and/or {self.exptime_key} keywords from FITS header')
            tmid = None
        return tmid

    @property
    def pgcenter(self):
        return PGCenter(self.obj_center,
                        self.desired_center,
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

class PGDNoCenter(PGData):
    """Sets :prop:`obj_center` to invalid value"""
    pass

class PGDCentered(PGData):
    """Sets :prop:`obj_center` to :prop:`desired_center`"""
    @property
    def obj_center(self):
        self._obj_center = self.desired_center
        self.quality = 10
        return self._obj_center

class PGDOffset(PGData):
    """Offsets :prop:`obj_center` by :prop:`center_offset`"""

    def __init__(self, *args,
                 center_offset=None,
                 **kwargs):
        super().__init__(*args, **kwargs)
        self._center_offset = None
        self.center_offset = center_offset

    @classmethod
    def read(cls, *args,
             center_offset=None,
             **kwargs):
        pgd = super(PGDOffset, cls).read(*args, **kwargs)
        pgd.__init__(from_read=True,
                     center_offset=center_offset)
        return pgd

    @property
    def center_offset(self):
        return self._center_offset

    @center_offset.setter
    def center_offset(self, value):
        print('in center_offset.setter')
        old_center_offset = self.center_offset
        if value is None:
            self.center_offset = (0,0)
        else:
            self._center_offset = np.asarray(value)
        if np.any(self._center_offset != old_center_offset):
            # Prepare to recalculate self._obj_center if we change the offset
            self._obj_center = None
        print(self.center_offset)

    @property
    def obj_center(self):
        if self._obj_center is not None:
            return self._obj_center
        self._obj_center = super().obj_center + self.center_offset
        return self._obj_center


class MaxImPGData(PGData):
    """MaxIM DL adjustments to the PGData class"""

    @property
    def binning(self):
        """Image binning in Y,X order"""
        binning = np.asarray((self.meta['YBINNING'],
                              self.meta['XBINNING']))
        return np.asarray(binning)
        
    @property
    def subframe_origin(self):
        """Subframe origin in *unbinned* pixels with full CCD origin = (0,0).  Y,X order"""
        subframe_origin = np.asarray((self.meta['YORGSUBF'],
                                      self.meta['XORGSUBF']))
        subframe_origin *= self.binning
        return subframe_origin

#pgc = PGCenter()
#pgc = PGCenter((1,2), (3,4))
#pgd = PGData()
fname = '/data/io/IoIO/raw/2020-07-15/HD87696-0016_Na_off.fit'

ccd = CCDData.read(fname, unit='adu')

pgd = PGData.read(fname)

rpgd = PGData(pgd.data, meta=pgd.meta)

dspgd = PGData.read(fname, desired_center=(100,100))
print(dspgd.obj_center, dspgd.desired_center)

dspgd.write('/tmp/test.fits', overwrite=True)

pgdc = PGDCentered.read(fname)
print(pgdc.obj_center, pgdc.desired_center)

pgdo = PGDOffset(pgd.data, meta=pgd.meta, center_offset=(20,10))

#pgdo = PGDOffset.read(fname)
#print(pgdo.obj_center, pgdo.desired_center)
#
#pgdo = PGDOffset.read(fname, center_offset=(20,10))
#print(pgdo.obj_center, pgdo.desired_center)

print('done')
