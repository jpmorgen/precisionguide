"""Define the base data containing classes for precisionguide system

These are designed to be portable across any platform, since they do
not depend on the specifics of the telescsope control system.  Thus,
data can be recorded using this code on, e.g., a Windows system and
uploaded to any other platform for detailed analysis using this same
code.

"""
import inspect 
import numpy as np
from scipy.signal import convolve2d
from scipy.ndimage.measurements import center_of_mass

from astropy.nddata import CCDData
from astropy import units as u
from astropy.time import Time
from astropy.utils import lazyproperty

_NotFound = object()

class pglazyproperty(lazyproperty):
    """Caching property decorator with auto setting and None reset

    Property in a Python object is usually passive:

    >>> class MyClass
    >>> self.myprop = 'a'
    >>> b = myfunct(self.myprop)

    The property decorator enables property to be active, such that a
    method return value can be considered as property.  The function
    that implements this is the "getter," with the syntax simplified:

    >>> @property
    >>> def mycalcprop(self):
    >>>     a = 1 + 2
    >>>     return a
    >>> print(self.mycalcprop)
    3

    Note that any time ``self.mycalcpropc`` is used, its return value
    is substituted just as if that was the real property.
    Importantly, with the configuration above, Python will not let you
    set the property.  In other words, without another step,
    self.mycalcprop cannot go on the left side of an ``=``.  The
    @mycalcproc.setter decorator provides the machinery for that
    next step.

    >>> @mycalcproc.setter
    >>> def mycalcproc(self, value):
    >>>     return value
    (no par is great if
    you just want a "one way" property

@property decorator refinement that enables calculations to be done
    once and stored without a setter (`astropy.utils.pglazyproperty`
    contribution) and enable the ``@<property>.setter`` to easily
    modify the user's input value for storage in the caching system
    (local contribution).  A further local refinement is to allow both
    initial and subsequent assignment of values to None that leave the
    system in a state that will run the getter's machinery the next
    time the property is queried.  If the getter takes a long time to
    calculate ``None,`` well, there will be no time savings, but at
    least the system has been arranged so there is no recursion error
    when ``None`` is the first value.  Example use: enable a setter to be
    given a ``tuple`` but covert it to a ``numpy.array`` for storage
    (see also :class:`pglazycoordproperty`).  Use ``del
    obj.<property>`` or set the property to None to delete the stored
    value of <property> and redo the calculation.

    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._none_proxy = 'None proxy'
    def __get__(self, obj, owner=None):
        try:
            obj_dict = obj.__dict__
            val = obj_dict.get(self._key, _NotFound)
            if val is self._none_proxy:
                return None
            if val is _NotFound or val is None:
                # I have no idea why self._lock not in the namespace
                #with self._lock:
                    # Check if another thread beat us to it.
                    #val = obj_dict.get(self._key, _NotFound)
                    #if val is _NotFound or val is None:
                        val = self.fget(obj)
                        if val is None:
                            obj_dict[self._key] = self._none_proxy
                            return None
                        obj_dict[self._key] = val                
            return val
        
        except AttributeError:
            if obj is None:
                return self
            raise

    def __set__(self, obj, val):
        print('setter')
        obj_dict = obj.__dict__
        if self.fset:
            ret = self.fset(obj, val)
            if ret is not None:
                obj_dict[self._key] = ret
                return
        obj_dict[self._key] = val        

class pglazycoordproperty(pglazyproperty):
    """Stores property values as np.arrays (e.g., converts tuples)
    Also enables recalculation when the value is set to None (or
    value could be deleted as in lazyproperty"""

    def __set__(self, obj, val):
        obj_dict = obj.__dict__
        if val is None:
            if obj_dict.get(self._key) is not None:
                # Don't delete value in __dict__ unless we have one!
                del obj_dict[self._key]
            return
        val = np.asarray(val)
        if val.shape != (2,):
            raise ValueError(f'value "{val}" is not a sensible 2D coordinate')
        super().__set__(obj, val)

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
    @pglazycoordproperty
    def obj_center(self):
        pass
        
    @pglazycoordproperty
    def desired_center(self):
        pass

    @pglazyproperty
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
        
    @pglazycoordproperty
    def tmid(self):
        pass


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
                 obj_center=None,
                 desired_center=None,
                 quality=0,
                 date_obs_key='DATE-OBS',
                 exptime_key='EXPTIME',
                 darktime_key='DARKTIME',
                 **kwargs):
        # Order is important!  Handle superclass first then subclass
        # for the reasons detailed below.
        #
        # Structure __init__ so that both standard instantiation and
        # .read classmethod can use the same code.  This enables us to
        # have one place where we separate out our subclass keywords
        # from the superclass keywords (**kwargs)

        # Handle super first to ensure image and meta property are
        # installed into the object before we invoke any code later in
        # the __init__ that might assume its presence
        if unit is None:
            unit = raw_unit
        if from_read:
            # Instantiate a new object with the super().read
            # classmethod and copy its property into our current
            # object.  Copying property over rather than just
            # reassigning "self" is necessary because when an object
            # is instantiated, an object reference is created that is
            # permanently associated with the object.  Simple
            # reassignment of self (e.g. self = ccd) does not have the
            # expected effect -- it just changes the place self points
            # to in the routine the assignment was made (e.g.,
            # __init__ in this case) and then when the routine exits,
            # the original "self" is used.  Note also that when we
            # call our super().read method and the object is
            # instantiated, this __init__ is read all over again,
            # since a new object is being created.  Since we don't
            # pass any parameters relevant to the subclass when we do
            # our super().read, all subclass property remains set to
            # the defaults.  When we __dict__.update our local
            # property, we copy over this freshly initialized property
            # to our current property.  Below, we re-initialize the
            # property with the actual keywords we provided to the
            # original .read call
            ccd = super().read(*args, unit=unit, **kwargs)
            self.__dict__.update(ccd.__dict__)
        else:
            super().__init__(*args, unit=unit, **kwargs)
        ## Generally the object will calculate these
        #if obj_center is not None:
        #    self.obj_center = obj_center
        #if desired_center is not None:
        #    self.desired_center = desired_center
        #if quality is not None:
        #    self.quality = quality
        # Generally the object will calculate these
        self.obj_center = obj_center
        self.desired_center = desired_center
        self.quality = quality
        self.date_obs_key = date_obs_key
        self.exptime_key = exptime_key
        self.darktime_key = darktime_key

    @classmethod
    def read(cls, *args, **kwargs):
        return cls(*args, from_read=True, **kwargs)

    @pglazycoordproperty
    def obj_center(self):
        # Here is our fancy calculation that takes a long time.
        # DESIGN NOTE: obj_center and desired_center use 
        # the pglazyproperty system so that they can be set to an
        # arbitrary value at any time.  Recalculation is triggered by 
        # and stored as an np.array without recalculation.  In
        # other words, we calculate just once but we enable 
        obj_center = np.asarray((-99, -99))
        self.quality = 0
        return obj_center

    @pglazycoordproperty
    def desired_center(self):
        # Here is where the complicated desired_center calculation is
        # done:
        return np.asarray(self.data.shape)/2

    @pglazyproperty
    def quality(self):
        self.obj_center

    @pglazyproperty
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

    @pglazyproperty
    def data_unbinned(self):
        """Returns an unbinned version of data
        """
        return self.im_unbinned(self.data)

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

class NoCenterPGD(PGData):
    """Sets :prop:`obj_center` to invalid value and :prop:`quality` to 0`"""
    pass

class CenteredPGD(PGData):
    """Sets :prop:`obj_center` to :prop:`desired_center`"""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._old_desired_center = None
    

    @pglazycoordproperty
    def obj_center(self):
        self.quality = 10
        return self.desired_center

    @pglazycoordproperty
    def desired_center(self):
        return super().desired_center

    @desired_center.setter
    def desired_center(self, value):
        # Reset obj_center for proper recalculation
        self.obj_center = None
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
        # pglazycoordproperty such that we want to delete that
        # property for recalculation in the event center_offset
        # changes.  Particularly considering we end up in np.arrays
        # and the += operator increments the object numpy object
        # itself (e.g. pgd.center_offset += (10,10)), this is just a
        # little too complex for the pglazycoordproperty system to
        # handle, so do it the olde fashionede waye

        # This copy is key for avoiding subtle assign by reference bugs
        if center_offset is None:
            center_offset = (0,0)
        self._old_center_offset = center_offset
        self.center_offset = center_offset
        #self._center_offset = np.asarray(center_offset)

    @pglazycoordproperty
    def calculated_center(self):
        # It is important this be its own name, since the lazyproperty
        # system works on the __dict__ key in the current object.  So
        # it doesn't "stick" from the super, but it does here under a
        # new name
        return super().obj_center

    #@property
    #def center_offset(self):
    #    return self._center_offset

    @pglazycoordproperty
    def _old_center_offset(self):
        return np.asarray((0,0))

    @pglazycoordproperty
    def center_offset(self):
        return np.asarray((0,0))

    @center_offset.setter
    def center_offset(self, value):
        print(value)
        if value is None:
            value = (0,0)
        value = np.asarray(value)
        print(value)
        if not np.all(self._old_center_offset == value):
            self._old_center_offset = value.copy()
            self.obj_center = None
        return value

    @pglazycoordproperty
    def obj_center(self):
        return self.calculated_center + self.center_offset

class MaxPGD(PGData):
    @pglazyproperty
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

    @pglazyproperty
    def obj_center(self):
        # Still working on this
        return super().obj_center + self.center_offset
    
class CenterOfMassPGD(PGData):
    """Use simple center-of-mass calculation to find center of
    brightest object.  Works best if pixels not likely to be part of
    the source are set to zero."""

    @pglazyproperty
    def obj_center(self):
        # Still working on this
        return super().obj_center + self.center_offset

class MaxImPGD(PGData):
    """Handles MaxIM DL :prop:`binning` and :prop:`subframe_origin`"""

    # It is important to just read these once, since we may modify
    # them in the object before write
    @pglazyproperty
    def binning(self):
        """Image binning in Y,X order"""
        binning = np.asarray((self.meta['YBINNING'],
                              self.meta['XBINNING']))
        return np.asarray(binning)
        
    @binning.setter
    def binning(self, value=None):
        return np.asarray(value)

    @pglazyproperty
    def subframe_origin(self):
        """Subframe origin in *unbinned* pixels with full CCD origin = (0,0).  Y,X order"""
        subframe_origin = np.asarray((self.meta['YORGSUBF'],
                                      self.meta['XORGSUBF']))
        subframe_origin *= self.binning
        return subframe_origin

    @subframe_origin.setter
    def subframe_origin(self, value=None):
        return np.asarray(value)

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
#print(pgc.obj_center)
#pgc.obj_center = (1,1)
#pgc.obj_center = 1
#print(pgc.obj_center)
#pgc = PGCenter((1,2), (3,4))
#pgd = PGData()
fname = '/data/io/IoIO/raw/2020-07-15/HD87696-0016_Na_off.fit'
#
class MyPGD(OffsetPGD, CenteredPGD, PGData):
    pass

#pgd = PGData.read(fname)
#print(pgd.obj_center, pgd.desired_center)
#pgd.obj_center = (1,1)
#print(pgd.obj_center, pgd.desired_center)
#pgd = PGData.read(fname, obj_center=(2,2))
#print(pgd.obj_center, pgd.desired_center)
#pgd.desired_center = (1,1)
#del pgd.desired_center
#pgd.obj_center = None
#print(pgd.obj_center, pgd.desired_center)
pgd = MyPGD.read(fname)
print(pgd.obj_center, pgd.desired_center, pgd.center_offset)
pgd.center_offset += np.asarray((10, -10))
pgd.center_offset += (10, -10)
print(pgd.obj_center, pgd.desired_center, pgd.center_offset)
pgd.center_offset = None
pgd.desired_center = (3,3)
#pgd.center_offset += (10, -30)
#pgd.center_offset = (10, -10)
print(pgd.obj_center, pgd.desired_center, pgd.center_offset)
#pgd.obj_center = (0,0)
#print(pgd.obj_center, pgd.desired_center, pgd.center_offset)
#pgd.obj_center = None
#del pgd.obj_center
#pgd.calculated_center = (0,0)
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
#class MyPGD(OffsetPGD, MaxPGD):
#    pass
#pgd = MyPGD.read(fname, center_offset=(-10,20))
#print(pgd.obj_center)
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

