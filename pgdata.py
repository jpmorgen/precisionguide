"""Define the base data containing classes for precisionguide system

These are designed to be portable across any platform, since they do
not depend on the specifics of the telescsope control system.  Thus,
data can be recorded using this code on, e.g., a Windows system and
uploaded to any other platform for detailed analysis using this same
code.

"""

from copy import deepcopy

import numpy as np

from scipy.signal import convolve2d
from scipy.ndimage.measurements import center_of_mass

from astropy import log
from astropy.io import fits
from astropy import units as u
from astropy.time import Time

from ccdmultipipe.utils.ccddata import FbuCCDData

INVALID_CENTER = (-99, -99)
INVALID_CENTER_QUALITY = 0

# These are generic keywords that need to get tailored
DATE_OBS_KEY = 'DATE-OBS'
EXPTIME_KEY = 'EXPTIME'
EXPTIME_UNIT = u.s
DARKTIME_KEY = 'DARKTIME'

_NotFound = object()

###########################
# Decorators              #
###########################
class pgproperty(property):
    """Caching property decorator with auto setting and None reset

    . myprop = None resets system and forces the getter to run

    . when getter runs and returns a non-None value, that value is
    cached and returned instead of running the getter

    . myprop = 'some_value' overrides the getter -- 'some_value' is
    cached and returned on subsequent queries of myprop

    . No setter is required to set the property, but if one is
    provided, its return value is used as the cached value that
    overrides the getter (i.e., no knowledge of the internal workings
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
        . tuple: like 0 but shape of value must match shape_check or
        a ValueError is raised

    Inspired by `astropy.utils.lazyproperty` and `property_manager
    <https://github.com/xolox/python-property-manager>`

    .NOTE: This cache operates on the top-level property.  In a
    complicated multi-inheritance system with inter-connected property
    like `pgdata`, some time savings in calculation of quantities
    could be achieved by caching property at each inheritance level
    individually, say by making the `key` a MRO level-unique string,
    e.g. property name, `self.__class__.__name__` and module name.
    Then a None reset would only reset the levels from the top down to
    the MRO at which the property was set to None.  This guarantees
    the children get a freshly calculated value of the affected
    properties and acknowledges that the parents don't care that the
    change was made, since they never depended on the relationship
    between the property in the first place.  Or if they did, they
    should be super()ed into the calculation and the most senior
    setter concerned about the interrelationship would also be calling
    for the reset of the affected other property.  The side-affect of
    caching all of the intermediate results would be more memory use,
    since all of the intermediate property would have long-lived
    references.  For the pgdata coordiate-oriented stuff, this is not
    a big deal, but it is not generally appropriate.

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
        if self.shape_check is None or val is None:
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
        val = self.npval(val)
        obj_dict[self._key] = val

    def __delete__(self, obj):
        if self.fdel:
            self.fdel(obj)
        obj.__dict__.pop(self._key, None)    # Delete if present


class pgcoordproperty(pgproperty):
    shape_check=(2,)


#######################
# Primary Objects     #
#######################

def center_quality_checker(value):
    if value is None:
        return
    if not isinstance(value, int) or value < 0 or value > 10:
        raise ValueError('quality must be an integer value from 0 to 10')
    return value

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

    center_quality : int
        0 -- 10 ranking of quality of ``obj_center``
        determination, with 1 = very bad, 10 = very good, 0 = invalid

    tavg : `~astropy.time.Time`
        Average time of observation.  See `:prop:~precisionguide.PGData.tavg`

    """
    def __init__(self,
                 obj_center=None,
                 desired_center=None,
                 center_quality=None,
                 tavg=None):
        self.obj_center = obj_center
        self.desired_center = desired_center
        self.center_quality = center_quality
        self.tavg = tavg

    # Our decorator does everything we need :-)
    @pgcoordproperty
    def obj_center(self):
        pass
        
    @pgcoordproperty
    def desired_center(self):
        pass

    @pgproperty
    def center_quality(self):
        pass
        
    @center_quality.setter
    def center_quality(self, value):
        """Unset center_quality translates to 0 center_quality"""
        self._center_quality = center_quality_checker(value)
    
    @pgproperty
    def tavg(self):
        pass

# We want to use astropy's CCDData to store our 


# Pretty much all CCD-like detectors present their raw data in adu.
# In the event this needed to change for a particular system, insert a
# subclass redefining raw_unit anywhere on the PGData inheritance
# chain.  That raw_unit will override this one.  Note, doing this as a
# class variable is necessary because of the CCDData read classmethod
class ADU():
    """Simple class to define the class variable `fallback_unit` =
    `~astropy.units.adu` for classes using the
    :class:`ccdmultipipe.FbuCCDData` system

    """
    fallback_unit = u.adu

class PGData(ADU, FbuCCDData):
    """Base class for image data in the `precisionguide` system

    This class stores CCD data and defines methods and property to
    calculate/store four primary quantities: :prop:`obj_center`,
    :prop:`desired_center`, :prop:`center_quality`, and :prop:`tavg`, as
    described in the documentation.

    CCD data are stored using `astropy.nddata.CCDData` as its base
    class with the addition of a class,
    :class:`ccdmultipipe.FbuCCDData`, that ensures the
    `~astropy.nddata.CCDData` will always be defined with a
    :class:`astropy.units.Unit`.  Since all or nearly all CCD-like
    detectors present their raw data in `~astropy.units.adu`, this is
    the default :class:`~astropy.units.Unit` assumed by
    :class:`PGData`.  In the event this needs to change for a
    particular system, a class similar to :class:`ADU` (e.g. "Count")
    could be defined and inserted into the inheritance chain


    superclass and adds properties and methods that calculate/store
    four quantities: :prop:`obj_center`, :prop:`desired_center`,
    :prop:`center_quality`, and :prop:`tavg`.  These four quantities are
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
                 data,
                 obj_center=None,
                 desired_center=None,
                 center_quality=None,
                 recalculate=False,
                 subframe_origin=None,
                 binning=None,
                 date_obs_key=DATE_OBS_KEY,
                 exptime_key=EXPTIME_KEY,
                 exptime_unit=EXPTIME_UNIT,
                 darktime_key=DARKTIME_KEY,
                 obs_location=None,
                 copy=False,
                 **kwargs):
        # Pattern after NDData init but skip all the tests
        if isinstance(data, PGData):
            # Sigh.  We have to undo the convenience of our pgproperty
            # lest we trigger the calculation of property, which leads
            # to recursion problems
            obj_dict = data.__dict__
            obj_center = obj_dict.get('obj_center')
            desired_center = obj_dict.get('desired_center')
            center_quality = data._center_quality
            subframe_origin = data.subframe_origin
            binning = data.binning
            date_obs_key = data.date_obs_key
            exptime_key = data.exptime_key
            exptime_unit = data.exptime_unit
            darktime_key = data.darktime_key
            obs_location = obs_location
        if copy:
            obj_center = deepcopy(obj_center)
            desired_center = deepcopy(desired_center)
            center_quality = deepcopy(center_quality)
            subframe_origin = deepcopy(subframe_origin)
            binning = deepcopy(binning)
            date_obs_key = deepcopy(date_obs_key)
            exptime_key = deepcopy(exptime_key)
            exptime_unit = deepcopy(exptime_unit)
            darktime_key = deepcopy(darktime_key)            
            obs_location = deepcopy(obs_location)            

        super().__init__(data, copy=copy, **kwargs)
        self.recalculate = recalculate # May be not appropriate at this level
        self.obj_center = obj_center
        self.desired_center = desired_center
        self._center_quality = center_quality
        self.subframe_origin = subframe_origin
        self.binning = binning
        self.date_obs_key = date_obs_key
        self.exptime_key = exptime_key
        self.exptime_unit = exptime_unit
        self.darktime_key = darktime_key
        self.obs_location = obs_location

    def _init_args_copy(self, kwargs):
        """The NDData _slice and _arithmetic methods create new class instances by running the __init__ method with the appropriate kwarg to copy over relevant property.  Thus we need to copy all our property into a keyword dictionary too"""
        # Grab calculated property directly from our pgproperty
        # decorator system so it doesn't trigger calculations!
        obj_dict = self.__dict__
        obj_center = obj_dict.get('obj_center')
        desired_center = obj_dict.get('desired_center')
        kwargs['obj_center'] = obj_center
        kwargs['desired_center'] = desired_center
        kwargs['center_quality'] = self._center_quality
        kwargs['subframe_origin'] = self.subframe_origin
        kwargs['binning'] = self.binning
        kwargs['date_obs_key'] = self.date_obs_key
        kwargs['exptime_key'] = self.exptime_key
        kwargs['darktime_key'] = self.darktime_key
        kwargs['obs_location'] = self.obs_location
        return kwargs
        
    def _arithmetic(self, *args, **kwds):
        """Insert our property when new class instances are created using arithmetic"""
        result, kwargs = super()._arithmetic(*args, **kwds)
        kwargs = self._init_args_copy(kwargs)
        return result, kwargs

    def _slice(self, item):
        """Override NDSlicingMixin definition to move subframe origin"""
        kwargs = super()._slice(item)
        kwargs = self._init_args_copy(kwargs)
        yslice = item[0]
        xslice = item[1]
        yorg = yslice.start or 0
        xorg = xslice.start or 0
        kwargs['subframe_origin'] = (yorg, xorg)
        return kwargs

    @pgcoordproperty
    def obj_center(self):
        """Center of object, (Y, X) in pixel

        Coordinates are referenced to the image stored in the
    :prop:`data` property.  This image may be a binned subframe of
    full detector.  Use :meth:`unbinned(self.obj_center)` to obtain
    the coordinates in raw detector pixels.  Quality of center
    determination must be set as well.  Base class uses out-of-bounds
    value of `precisionguide.INVALID_CENTER` for center and
    `precisionguide.INVALID_CENTER_QUALITY` for center_quality

        Results are stored using the :class:`pgcoordproperty`
        decorator system.  See documentation of that class for
        explanation of features.

        """
        obj_center = INVALID_CENTER
        self.center_quality = INVALID_CENTER_QUALITY
        return obj_center

    @pgcoordproperty
    def desired_center(self):
        """Desired center of object (Y, X).  Base class uses center of
    data array.  As with :prop:`obj_center`, :prop:`desired_center` is
    referenced to the data stored in :prop:`data`

        Results are stored using the :class:`pgcoordproperty`
        decorator system.  See documentation of that class for
        explanation of features.

        """
        # Here is where the complicated desired_center calculation is
        # done:
        npshape = np.asarray(self.data.shape)
        desired_center = npshape/2
        return desired_center

    @property
    def center_quality(self):
        """Quality of center determination on a 0 to 10 integer scale.

  Quality should always be set in the :prop:`obj_center` setter

        """
        if self._center_quality is None:
            self.obj_center
        return self._center_quality

    @center_quality.setter
    def center_quality(self, value):
        self._center_quality = center_quality_checker(value)

    @pgproperty
    def obs_location(self):
        """Return `~astropy.coordinates.EarthLocation` of observatory.  NOTE:
        for the basic precisionguide system to work, this property is
        not necessary because only time deltas between observations at
        the observatory are used in calculations.  Setting this
        property in a subclass is desirable if
        `~:prop:precisionguide.PGData.tavg` of that subclass class is
        used in calculations involving light travel time between the
        observatory and a distant point (e.g. the solar system
        barycenter).

        """
        return None
        # SAMPLE CODE:
        #from astropy.coordinates import EarthLocation
        #lon = hdr.get('long-obs')
        #lat = hdr.get('lat-obs')
        ## Just use an average elevation if none is available, since that
        ## is probably better than default of elevation = 0
        ## https://www.quora.com/What-is-the-average-elevation-of-Earth-above-the-ocean-including-land-area-below-sea-level-What-is-the-atmospheric-pressure-at-that-elevation
        #alt = hdr.get('alt-obs') or 800 * u.m
        #return EarthLocation.from_geodetic(lon, lat, alt)

    @pgproperty
    def tavg(self):
        """`~astropy.time.Time` Average time of observation.  For simple
        exposures, this is the midpoint.  For exposures integrated
        from several segments one might imagine this would be the
        average of the exposure efficiency as a function of time.
        However Rots et al (A&A 574, A36 2015) note no standard is
        defined.  This can turn into the DATE-AVG FITS keyword,
        although astropy.wcsparm deletes that when a
        `~astropy.nddata.CCDData` is read in.  The astropy standard is
        to assume the DATE-OBS keyword represents the start of the
        exposure.  There is no standard for exposure time, hence the
        `:prop:precisionguide.PGData.exptime` property.  For timing
        calculations involving two distant points and/or moving
        frames, the observatory location should be specified via
        `:prop:precisionguide.PGData.obs_location`

        """ 
        timesys = self.meta.get('TIMESYS') or 'UTC'
        tavg_str = self.meta.get('DATE-AVG')
        if tavg_str is not None:
            return Time(tavg_str, format='fits',
                        scale=timesys.lower(),
                        location=self.obs_location)
        try:
            exptime = self.meta[self.exptime_key.lower()]
            exptime *= self.exptime_unit
            dateobs_str = self.meta[self.date_obs_key.lower()]
            tbeg = Time(dateobs_str, format='fits',
                        scale=timesys.lower(),
                        location=self.obs_location)
            return tbeg + exptime/2
        except:
            log.error(f'Cannot read sufficient combination of '
                      '{self.date_obs_key}, {self.darktime_key} '
                      'and/or {self.exptime_key} keywords from FITS '
                      'header to establish tavg')
            return None

    @tavg.setter
    def tavg(self, val):
        if not isinstance(val, Time):
            raise ValueError('tavg must be of type `astropy.time.Time`')

    @pgcoordproperty
    def binning(self):
        """Image binning in Y,X order.

         NOTE: this needs to be overridden in a subclass with the
         actual FITS binning keywords used (which are unfortunately
         not in any FITS standard).  E.g.:
         >>> binning = (self.meta['YBINNING'],
         >>>            self.meta['XBINNING'])
         >>> return binning
        """
        pass

        #binning = (1,1)
        #return binning
        
    @pgcoordproperty
    def subframe_origin(self):
        """Subframe origin in *unbinned* pixels with full CCD origin =
    (0,0).  Y,X order
         NOTE: this needs to be overridden in a subclass with the
         actual FITS binning keywords used (which are unfortunately
         not in any FITS standard).  E.g.:
         >>> subframe_origin = (self.meta['YORGSUBF'],
         >>>                    self.meta['XORGSUBF'])
         >>> subframe_origin = np.asarray(subframe_origin)
         >>> subframe_origin *= self.binning
         >>> return subframe_origin
        """
        pass
        #subframe_origin = (0,0)
        #return subframe_origin

    def coord_unbinned(self, coords):
        """Returns pixel coords referenced to full CCD given internally stored binning/subim info"""
        coords = np.asarray(coords)
        unbinned = self.binning * coords + self.subframe_origin
        return unbinned.astype(int)

    def x_unbinned(self, coords):
        """Returns x coord referenced to full CCD given internally stored binning/subim info"""
        coords = np.asarray(coords)
        unbinned = self.binning[1] * coords + self.subframe_origin[1]
        return unbinned.astype(int)
        
    def y_unbinned(self, coords):
        """Returns y coord referenced to full CCD given internally stored binning/subim info"""
        coords = np.asarray(coords)
        unbinned = self.binning[0] * coords + self.subframe_origin[0]
        return unbinned.astype(int)

    def coord_binned(self, coords, limit_edges=False):
        """Assuming coords are referenced to full CCD, return location in binned coordinates relative to the subframe origin.  Limit_edges=True sets output value(s) that would otherwise fall outside the binned image to the appropriate boundary value """
        coords = np.asarray(coords)
        binned_coords = (coords - self.subframe_origin) / self.binning
        if not limit_edges:
            return binned_coords
        s = np.asarray(self.shape)
        if len(coords.shape) == 1:
            binned_coords[0] = np.min((s[0], binned_coords[0]))
            binned_coords[1] = np.min((s[1], binned_coords[1]))
            return binned_coords.astype(int)
        lidx = np.flatnonzero(binned_coords[0, :] > s[0])
        binned_coords[0, lidx] = s[0]
        lidx = np.flatnonzero(binned_coords[1, :] > s[1])
        binned_coords[1, lidx] = s[1]
        return binned_coords.astype(int)

    def x_binned(self, xs, **kwargs):
        """Assuming x coords are referenced to full CCD, return location in binned coordinates relative to the subframe origin.  Limit_edges=True sets output value(s) that would otherwise fall outside the binned image to the appropriate boundary value """
        if np.isscalar(xs):
            unbinned_coords = np.asarray((-99, xs))
            binned_coords = self.coord_binned(unbinned_coords, **kwargs)
            xs = binned_coords[1]
        else:
            xs = np.asarray(xs)
            ys = len(xs)*(-99,)
            unbinned_coords = list(zip(ys, xs))
            binned_coords = self.coord_binned(unbinned_coords, **kwargs)
            xs = binned_coords[:, 1]
        return xs
       
    def y_binned(self, ys, **kwargs):
        """Assuming x coords are referenced to full CCD, return location in binned coordinates relative to the subframe origin.  Limit_edges=True sets output value(s) that would otherwise fall outside the binned image to the appropriate boundary value """ 
        if np.isscalar(ys):
            unbinned_coords = np.asarray((ys, -99))
            binned_coords = self.coord_binned(unbinned_coords, **kwargs)
            ys = binned_coords[0]
        else:
            ys = np.asarray(ys)
            xs = len(ys)*(-99,)
            unbinned_coords = list(zip(ys, xs))
            binned_coords = self.coord_binned(unbinned_coords, **kwargs)
            ys = binned_coords[:, 0]
        return ys
       
    def im_unbinned(self, a):
        """Returns an unbinned version of a.  a must be same shape as self.data
        NOTE: to preserve flux, divide by (np.prod(self.binning)
        """
        if a is None:
            return None
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
    def self_unbinned(self):
        """Returns an unbinned version of this object.  NOTE: just as with the
    primary object, if this object needs to be modified, it should be
    copied first, particularly since, in the case the image is not
    binned, it might just point to self

        """
        if (self.binning.sum() == 2
            and self.subframe_origin.sum() == 0):
            # We are already unbinned, don't bother to copy
            return self
        # If we made it here, we need to unbin
        self_unbinned = deepcopy(self)
        self_unbinned.data = self.im_unbinned(self.data)
        self_unbinned.mask = self.im_unbinned(self.mask)
        self_unbinned.uncertainty = self.im_unbinned(self.uncertainty)
        self_unbinned.binning = np.asarray((1,1))
        self_unbinned.subframe_origin = np.asarray((0,0))
        # --> NOTE we are not messing with any other metadata
        return self_unbinned

    @property
    def pgcenter(self):
        """Returns a :class:`PGCenter` object with the *unbinned* obj and
    desired center coordinates.  Working in *unbinned* coordinates is
    essential for the internal workings of the precisionguide
    telescope control system, but rather confusing and abstract for
    the user in any other context.

        """
        # I may want to also export the binning information for
        # display purposes in precisionguide
        return PGCenter(self.coord_unbinned(self.obj_center),
                        self.coord_unbinned(self.desired_center),
                        self.center_quality,
                        self.tavg)

    def _card_write(self):
        """Writes FITS header cards for obj and desired center coordinates in
        *binned* coordinates (i.e., as they would likely presesnt in
        image analysis software showing the image)

        """
        if not np.all(np.isfinite(self.obj_center)):
            self.obj_center = INVALID_CENTER
            log.warning(f'non-finite obj_center found, converting to '
                        f'{self.obj_center}')
        # Note pythonic y, x coordinate ordering
        self.meta['OBJ_CR0'] = (self.obj_center[1], 'Calculated object '
                                'center X')
        self.meta['OBJ_CR1'] = (self.obj_center[0], 'Calculated object '
                                'center Y')
        self.meta['DES_CR0'] = (self.desired_center[1], 'Desired center X')
        self.meta['DES_CR1'] = (self.desired_center[0], 'Desired center Y')
        self.meta['HIERARCH CENTER_QUALITY'] = (
            self.center_quality or INVALID_CENTER_QUALITY,
            'Quality on 0-10 scale '
            'of center determination')
        # As per https://www.aanda.org/articles/aa/pdf/2015/02/aa24653-14.pdf
        if not self.meta.get('DATE-AVG'):
            self.meta.insert('DATE-OBS',
                             ('DATE-AVG', self.tavg.value,
                              'midpoint of exposure, UT'),
                             after=True)

    def write(self, *args, **kwargs):
        self._card_write()
        super().write(*args, **kwargs)

class NoCenterPGD(PGData):
    """Sets :prop:`obj_center` to invalid value and :prop:`center_quality` to 0`"""
    pass

class FITSReaderPGD(PGData):
    """Read FITS keys to set defaults"""

    ## Inefficient way -- FITS header is read each time object is
    ## initialized even if it never uses the property
    #def __init__(self, *args,
    #             obj_center=None,
    #             desired_center=None,
    #             center_quality=0,
    #             recalculate=False,
    #             **kwargs):
    #    super().__init__(*args, **kwargs)
    #    try:
    #        assert not recalculate, ('Recalculation requested')
    #        log.debug('Trying to read center info from FITS header,')
    #        cx = self.meta['OBJ_CR0']
    #        cy = self.meta['OBJ_CR1']
    #        dx = self.meta['DES_CR0']
    #        dy = self.meta['DES_CR1']
    #        q = self.meta['CENTER_QUALITY']
    #        tavg_str = self.meta['DATE-AVG']
    #        self.obj_center = (cy, cx)
    #        self.desired_center = (dy, dx)
    #        self.center_quality = q
    #        self.tavg = Time(tavg_str, format='fits')
    #        log.info('Center info set from FITS header, use '
    #                 'recalculate=True to avoid')
    #    except Exception as e:
    #        log.debug(str(e))
    #        log.debug('Setting obj_center, desired_center, and center_quality '
    #                  'from object instantiation keywords')
    #        self.obj_center = obj_center
    #        self.desired_center = desired_center
    #        self.center_quality = center_quality

    @pgproperty
    def obj_center(self):
        """Center of object, (Y, X).  Quality of center determination must be
    set as well.  Base class uses out-of-bounds value (-99, -99) for
    center and 0 for center_quality

        Results are stored using the :class:`pgcoordproperty`
        decorator system.  See documentation of that class for
        explanation of features.

        """
        try:
            assert not self.recalculate, ('Recalculation requested')
            log.debug('Trying to read center info from FITS header,')
            cx = self.meta['OBJ_CR0']
            cy = self.meta['OBJ_CR1']
            obj_center = np.asarray((cy, cx))
        except Exception as e:
            log.debug(str(e))
            log.debug('Not setting center from FITS header,')
            obj_center = None
            self.center_quality = None
        return obj_center

    @pgproperty
    def desired_center(self):
        """Desired center of object (Y, X).  Base class uses center of
    data array.

        Results are stored using the :class:`pgcoordproperty`
        decorator system.  See documentation of that class for
        explanation of features.

        """
        try:
            assert not self.recalculate, ('Recalculation requested')
            log.debug('Trying to read desired center info from FITS header,')
            dx = self.meta['DES_CR0']
            dy = self.meta['DES_CR1']
            desired_center = np.asarray((dy, dx))
        except Exception as e:
            log.debug(str(e))
            log.debug('Not setting desired center from FITS header,')
            desired_center = None
        return desired_center

    @pgproperty
    def center_quality(self):
        """Quality of center determination.  center_quality should always be
        set in the obj_center setter

        """
        try:
            assert not self.recalculate, ('Recalculation requested')
            log.debug('Trying to read center_quality info from FITS header,')
            q = self.meta['CENTER_QUALITY']
            return int(q)
        except Exception as e:
            log.debug(str(e))
            log.debug('Not setting center_quality from FITS header,')
            # obj_center should set center_quality
            self.obj_center

    @center_quality.setter
    def center_quality(self, value):
        """Unset center_quality translates to 0 center_quality"""
        if value is not None:
            value = center_quality_checker(value)
        return value

#### In the end, I think CameraData is a bad idea, since this can all
#### go into and out of metadata.  Unit in principle could go there
#### too (e.g. BUNIT), however it is implemented as separate property
#### in the underlying CCDData/NDData.  No need to add that level of
#### complexity to these.  FITS card metadata property with comments
#### is good enough.
###
#### Although CameraData is intended to be for containing just camera
#### information, make it a subclass of FbuCCDData so that the raw_unit
#### can be properly inserted when CCD images are read in.
###
#### I envision either the camera file being use with a camera name and
#### method to initialize all of the propriety, or just a new class being
#### prepared that has the property hand-coded as class variables
###class CameraData(FbuCCDData):
###    raw_unit = u.adu
###    camera_data_file = None # some kind of file that would contain this info
###    camera_name = None
###    camera_description = None
###    full_naxis1 = None
###    full_naxis2 = None
###    satlevel = None
###    nonlinlevel = None
###    gain = None
###    readnoise = None
###
###    def __init__(self,
###                 *args,
###                 camera_data_file=None,
###                 camera_name=None,
###                 camera_description=None,
###                 raw_unit=None, # This will not affect the read classmethod
###                 full_naxis1=None,
###                 full_naxis2=None,
###                 satlevel=None,
###                 nonlinlevel=None,
###                 gain=None,
###                 readnoise=None,
###                 **kwargs):
###        self.camera_data_file = camera_data_file or self.camera_data_file
###        self.camera_name = camera_name or self.camera_name
###        self.camera_description = camera_description or self.camera_description
###        self.raw_unit = raw_unit or self.raw_unit
###        self.full_naxis1 = full_naxis1 or self.full_naxis1
###        self.full_naxis2 = full_naxis2 or self.full_naxis2
###        self.satlevel = satlevel or self.satlevel
###        self.nonlinlevel = nonlinlevel or self.nonlinlevel
###        self.gain = gain or self.gain
###        self.readnoise = readnoise or self.readnoise
###        super().__init__(*args, fallback_unit=self.raw_unit, **kwargs)
###
###    @classmethod
###    def read(cls, filename, 
###             raw_unit=None,
###             **kwargs):
###        """Use ``raw_unit`` instead of ``fallback_unit``."""
###        raw_unit = raw_unit or cls.raw_unit
###        return super(CameraData, cls).read(filename,
###                                       fallback_unit=raw_unit,
###                                       **kwargs)
###
###    def _card_write(self):
###        """Write FITS header cards for camera
###
###        """
###        self.meta['GAIN'] = (self.gain, f'CCD charge gain '
###                             '{self.gain.unit.to_str()}')
###        self.meta['SATLEVEL'] = (self.satlevel, f'CCD saturation level '
###                                 '{self.satlevel.unit.to_str()}')
###        self.meta['NONLIN'] = (self.nonlinlevel, f'CCD nonlinearity '
###                               'level {self.nonlinlevel.unit.to_str()}')
###        self.meta['RDNOISE'] = (self.readnoise, f'CCD readnoise '
###                               'level {self.readnoise.unit.to_str()}')
###
###class SX694(CameraData):
###    camera_name = 'SX694'
###    camera_description = 'Starlight Xpress Trius SX694 mono, 2017 model version'
###    raw_unit = u.adu
###    # naxis1 = fastest changing axis in FITS primary image = X in
###    # Cartesian thought
###    # naxis1 = next to fastest changing axis in FITS primary image = Y in
###    # Cartesian thought
###    full_naxis1 = 2750*u.pix
###    full_naxis2 = 2200*u.pix
###    # 16-bit A/D converter
###    satlevel = (2**16-1)*raw_unit
###    nonlinlevel = (42000 - 1811)*raw_unit
###    # Gain measured in /data/io/IoIO/observing/Exposure_Time_Calcs.xlsx.
###    # Value agrees well with Trius SX-694 advertised value (note, newer
###    # "PRO" model has a different gain value).  Stored in GAIN keyword
###    gain = 0.3 * u.electron/raw_unit
###    # Sample readnoise measured as per ioio.notebk
###    # Tue Jul 10 12:13:33 2018 MCT jpmorgen@byted 
###    # Readnoies is measured regularly as part of master bias creation and
###    # stored in the RDNOISE keyword.  This is used as a sanity check.
###    readnoise = 15.475665*raw_unit




class Ex(PGData):
    @pgcoordproperty
    def obj_center(self):
        # Pattern for all that want to read the FITS
        print('in obj_center')
        t = super().obj_center
        return t

class ExampleFITSReaderPGD(FITSReaderPGD):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        print('in init')              

    @pgcoordproperty
    def obj_center(self):
        # Pattern for all that want to read the FITS
        t = FITSReaderPGD(self).obj_center
        if  t is not None:
            return t
        # Real calculation starts here
        obj_center = self.desired_center
        self.center_quality = 10
        return obj_center

    @pgcoordproperty
    def desired_center(self):
        # Pattern for all that want to read the FITS
        t = FITSReaderPGD(self).desired_center
        if t is not None:
            return t
        # Real calculation starts here
        npshape = np.asarray(self.shape)
        desired_center = npshape/2
        return desired_center
    
    @pgproperty
    def center_quality(self):
        # Pattern for all that want to read the FITS
        t = FITSReaderPGD(self).center_quality
        if t is not None:
            return t
        self.obj_center
    

class CenteredPGD(PGData):
    """Sets :prop:`obj_center` to :prop:`desired_center`"""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._old_desired_center = None
    
    @pgcoordproperty
    def obj_center(self):
        self.center_quality = 10
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

    def _card_write(self):
        """Writes FITS header cards in *binned* coordinates (i.e., as they
would likely presesnt in image analysis software showing the image)

        """
        # Note pythonic y, x coordinate ordering
        self.meta['OFF_CR0'] = (self.center_offset[1], 'Offset from '
                                'orig obj center X')
        self.meta['OFF_CR1'] = (self.center_offset[0], 'Offset from '
                                'orig obj center Y')
        super()._card_write()

class MaxPGD(PGData):
    @pgcoordproperty
    def obj_center(self):
        self.center_quality = 6
        obj_center = np.unravel_index(np.argmax(self), self.shape)
        return obj_center

# --> Still working on this.  Looks like I am going to need a
# readnoise property in the object if I really want to do this.  Maybe
# we need a readnoise_keyword or something like that....  Very camera
# specific

# --> this could be property in pgdata and use biweight_local
class BackgroundPGD(PGData):
    """Hold background.  Might need this to expand to hold other
    things.  Util is a bit too general a name, but that is sort of
    what I would envision"""

    @pgproperty
    def background(self):
        return np.median(self)

class CenterOfMassPGD(BackgroundPGD):
    """Use simple center-of-mass calculation to find center of
    brightest object.  Works best if pixels not likely to be part of
    the source are set to zero."""

    @pgcoordproperty
    def obj_center(self):
        #bsub = self.subtract(self.background, handle_meta='first_found')
        #return center_of_mass(bsub.data)
        bsub = self.data - self.background
        self.center_quality = 6
        return center_of_mass(bsub)

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
    
class MaxImPGD(PGData):
    """Handles MaxIM DL :prop:`binning` and :prop:`subframe_origin`"""

    # It is important to just read these once, since we may modify
    # them in the object before write
    @pgcoordproperty
    def binning(self):
        """Image binning in Y,X order"""
        binning = (self.meta['YBINNING'],
                   self.meta['XBINNING'])
        binning = np.asarray(binning)
        return binning
        
    @pgcoordproperty
    def subframe_origin(self):
        """Subframe origin in *unbinned* pixels with full CCD origin = (0,0).  Y,X order"""
        subframe_origin = (self.meta['YORGSUBF'],
                           self.meta['XORGSUBF'])
        subframe_origin = np.asarray(subframe_origin)
        subframe_origin *= self.binning
        return subframe_origin

    def _card_write(self):
        """Write FITS card unique to MaxIMPGD"""
        # Note pythonic y, x coordinate ordering
        self.meta['XBINNING'] = self.binning[1]
        self.meta['YBINNING'] = self.binning[0]
        self.meta['XORGSUBF'] = self.subframe_origin[1]
        self.meta['YORGSUBF'] = self.subframe_origin[0]
        super()._card_write()

    def write(self, *args, **kwargs):
        self._card_write()
        super().write(*args, **kwargs)

if __name__ == "__main__":
    log.setLevel('DEBUG')
    #pgc = PGCenter()
    #pgc.obj_center = (1,1)
    #print(pgc.obj_center)
    ##pgc.obj_center = 1
    ###pgc = PGCenter((1,2), (3,4))
    #pgc = PGCenter([1,2], (3,4))
    #print(pgc.obj_center, pgc.desired_center)

    #pgd = PGData()
    fname = '/data/io/IoIO/raw/2020-07-15/HD87696-0016_Na_off.fit'
    
    class MyPGD(OffsetPGD, CenteredPGD, PGData):
        pass

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
    #
    #pgd = CenteredPGD.read(fname)
    #print(pgd.obj_center, pgd.desired_center)
    
    #pgd = MyPGD.read(fname)
    #print(pgd.obj_center, pgd.desired_center, pgd.center_offset)
    #pgd.center_offset += np.asarray((10, -10))
    #pgd.center_offset += (10, -10)
    #print(pgd.obj_center, pgd.desired_center, pgd.center_offset)
    #pgd.center_offset = None
    #pgd.desired_center = (3,3)
    #pgd.center_offset += (10, -30)
    #pgd.center_offset = (10, -10)
    #print(pgd.obj_center, pgd.desired_center, pgd.center_offset)
    #pgd.obj_center = (0,0)
    #print(pgd.obj_center, pgd.desired_center, pgd.center_offset)
    #pgd.obj_center = None
    ##del pgd.obj_center
    ###pgd.calculated_center = (0,0)
    #print(pgd.obj_center, pgd.desired_center, pgd.center_offset)
    #pgd.center_offset = pgd.center_offset + (10, -10)
    #print(pgd.obj_center, pgd.desired_center, pgd.center_offset)
    #pgd.center_offset += (10, -30)
    #print(pgd.obj_center, pgd.desired_center, pgd.center_offset)
    #pgd.write('/tmp/testfits.fits', overwrite=True)

    #class MyPGD(OffsetPGD, CenteredPGD, PGData):
    class MyPGD(OffsetPGD, CenteredPGD, MaxImPGD):#, SX694):
        pass
    #pgd = MyPGD.read(fname)
    #pgd.center_offset += (10, -30)
    #print(pgd.obj_center, pgd.desired_center, pgd.center_offset)
    #pgd.write('/tmp/testfits.fits', overwrite=True)

    #pgd = MyPGD(pgd)
    #print(pgd.obj_center, pgd.desired_center, pgd.center_offset)
    #pgd = MyPGD(pgd.data)
    #print(pgd.obj_center, pgd.desired_center, pgd.center_offset)

    #pgd = Ex.read(fname)
    #print(pgd.obj_center, pgd.desired_center)

    #pgd = FITSReaderPGD.read('/tmp/testfits.fits')
    #print(pgd.obj_center, pgd.desired_center)
    #
    #pgd = ExampleFITSReaderPGD.read('/tmp/testfits.fits')
    #print(pgd.obj_center, pgd.desired_center, pgd.center_quality)    
    #
    #pgd = ExampleFITSReaderPGD.read(fname)
    #print(pgd.obj_center, pgd.desired_center, pgd.center_quality)    

    #pgd = CenterOfMassPGD.read(fname)
    #print(pgd.obj_center, pgd.desired_center, pgd.center_quality)        
    #
    #off_filt = '/data/io/IoIO/raw/2018-01-28/R-band_off_ND_filter.fit'
    #pgd = CenterOfMassPGD.read(off_filt) 
    #print(pgd.obj_center, pgd.desired_center, pgd.center_quality)        

    #pgd.obj_center = (1,1)
    #pgd.desired_center = (1,1)
    #print(pgd.obj_center, pgd.desired_center)
    #del pgd.desired_center
    #print(pgd.obj_center, pgd.desired_center)
    #pgd.center_offset = (20. -10)
    #print(pgd.obj_center, pgd.desired_center, pgd.center_offset)

    #log.setLevel('DEBUG')
    #pgd = MaxPGD.read(fname)
    #print(pgd.obj_center)
    #wname = '/tmp/test_kwd_write.fits'
    #pgd.write(wname, overwrite=True)
    #pgd = MaxPGD.read(wname)
    #print(pgd.obj_center, pgd.desired_center)
    #pgd = PGData.read(wname, recalculate=True)
    #print(pgd.obj_center, pgd.desired_center)
    #class MyPGD(OffsetPGD, MaxPGD):
    #    pass
    ##pgd = MyPGD.read(fname, center_offset=(-10,20))
    ##print(pgd.obj_center, pgd.desired_center, pgd.center_offset)
    #
    #bname = '/data/io/IoIO/reduced/Calibration/2020-07-07_ccdT_-10.3_bias_combined.fits'
    #bpgd = MyPGD.read(bname)
    #epgd = MyPGD.read(bname, unit='adu')
    #ccd = FbuCCDData.read(bname)

    #pgd = MyPGD.read(fname)


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

    #bname = '/data/io/IoIO/reduced/Calibration/2020-07-07_ccdT_-10.3_bias_combined.fits'
    ##ccd = CCDData.read(bname)
    #
    #fname1 = '/data/Mercury/raw/2020-05-27/Mercury-0005_Na-on.fit'
    ##ccd = CCDData.read(fname1)
    #ccd = FbuCCDData.read(fname1, fallback_unit='adu')
    ##ccd = FbuCCDData.read(fname1, fallback_unit='aduu')
    #
    #ccd = FbuCCDData.read(fname1, unit='electron')

    fname1 = '/data/Mercury/raw/2020-05-27/Mercury-0005_Na-on.fit'
    pgd = CenterOfMassPGD.read(fname1)
    print(pgd.obj_center, pgd.desired_center)

#fname1 = '/data/IoIO/raw/20210310/HD 132052-S001-R001-C002-R.fts'
#pgd = MaxImPGD.read(fname1)
##print(pgd.meta)
#sub = pgd[0:100, 20:100]
#print(sub.subframe_origin)

#fname = '/data/IoIO/raw/2021-04_Astrometry/Main_Astrometry_East_of_Pier.fit'
#c = MaxImPGD.read(fname)
#print(c.binning)
#print(c.x_binned(6))
#print(c.x_binned(np.asarray((6, 12, 24, 48))))
#print(c.y_binned(6))
#print(c.y_binned(np.asarray((6, 12, 24, 48))))
