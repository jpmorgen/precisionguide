"""Define lightweight data containers for precisionguide system"""

import numpy as np

from astropy.io import fits

from astropy.time import Time

from astropy import units as u
from astropy.nddata import CCDData

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


class PGData():
    """Base class for image data in the `precisionguide` system

    This class contains an individual image data and calculates and/or
    stores four quantities: :prop:`obj_center`,
    :prop:`desired_center`, :prop:`quality`, and :prop:`tmid`.  These
    four quantities are intended to be returned in a :class:`PGCenter`
    object for subsequent lightweight storage and use in the
    precisionguide system.  Because precisionguide controls the
    absolute position of an object (or FOV center) on the CCD,
    :prop:`desired_center` and :prop:`obj_center` always read in
    *unbinned* pixel values referenced to the origin of the CCD
    itself.  Thus, the image input to :class:`PGData` must include
    both the image FITS header and image array(s)

    Parameters
    ----------
    input_im : `str`, `~astropy.io.fits.HDUList`, `~astropy.nddata.CCDData`, or `numpy.ndarray`

        Input FITS filename or object

    meta : `~astropy.io.fits.Header`
        FITS header.  Will replace any header read from file or
        object
        Default is ``None``
"""

    def __init__(self,
                 data,
                 meta=None,
                 desired_center=None,
                 date_obs_key='DATE-OBS',
                 exptime_key='EXPTIME',
                 darktime_key='DARKTIME'):
        self._data = data
        self._meta = meta
        self._data_unbinned = None
        self._desired_center = None
        self._obj_center = None
        self.desired_center = desired_center
        self.date_obs_key = date_obs_key
        self.exptime_key = exptime_key
        self.darktime_key = darktime_key

    @classmethod
    def read(cls, filename,
             meta=None,
             **kwargs):
        with fits.open(filename) as hdul:
            data = hdul[0].data
            meta = hdul[0].header
            return cls(data, meta=meta, **kwargs)



    #@classmethod
    #def read(cls, input_im, *args, **kwargs):
    #    return cls(input_im, *args, **kwargs)

    @property
    def data(self):
        return self._data
        
    @data.setter
    def data(self, value):
        # A little angst about allowing user to reset the data without
        # checking the center.  Generally, the object should not be
        # used as a container that gets reused for different images,
        # but it does seem convenient to be able to overwrite the 
        if self._data is not None and value.shape != self._data.shape:
            # Unfortunately, we can't check for subframe origin shift
            raise ValueError('New data array must be same shape as old array')
        self._data = value
        self._data_unbinned = None
        
    @property
    def meta(self):
        return self._meta
        
    @meta.setter
    def meta(self, value):
        self._meta = value
        
    @property
    def header(self):
        return self._meta
        
    @header.setter
    def header(self, value):
        self.meta = value

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

    def _card_write(self):
        # Note pythonic y, x coordinate ordering
        self.meta['DES_CR0'] = (self._desired_center[1], 'Desired center X')
        self.meta['DES_CR1'] = (self._desired_center[0], 'Desired center Y')
        self.meta['OBJ_CR0'] = (self._obj_center[1], 'Object center X')
        self.meta['OBJ_CR1'] = (self._obj_center[0], 'Object center Y')
        self.meta['QUALITY'] = (self.quality, 'Quality on 0-10 scale of center determination')
        self.meta['QUALITY'] = (self.quality, 'Quality on 0-10 scale of center determination')

    def write(self, *args, **kwargs):
        self._card_write()
        hdu = fits.PrimaryHDU(self.data, self.meta)
        hdu.writeto(*args, **kwargs)

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
        return PGCenter(self.obj_center, self.desired_center, self.quality)

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

class PGDOffsetCenter(PGDCentered):
    """Offsets :prop:`obj_center` by :prop:`center_offset`"""

    def __init__(self,
                 *args,
                 center_offset=None,
                 **kwargs):

        self._center_offset = None
        self.center_offset = center_offset
        super().__init__(*args, **kwargs)        

    @property
    def center_offset(self):
        return self._center_offset

    @center_offset.setter
    def center_offset(self, value):
        old_center_offset = self.center_offset
        if value is None:
            self.center_offset = (0,0)
        else:
            self._center_offset = np.asarray(value)
        if np.any(self._center_offset != old_center_offset):
            # Prepare to recalculate self._obj_center if we change the offset
            self._obj_center = None

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

class PGCCDData(PGData, CCDData):
    """Merges :class:`PGData and :class:`~astropy.nddata.CCDData`"""
    def __init__(self, *args,
                 meta=None,
                 unit=None,
                 raw_unit='adu',
                 desired_center=None,
                 date_obs_key='DATE-OBS',
                 exptime_key='EXPTIME',
                 darktime_key='DARKTIME',
                 **kwargs):
        if unit is None:
            unit = raw_unit
        CCDData.__init__(self, *args, meta=meta, unit=unit, **kwargs)
        PGData.__init__(self, self.data, meta=self.meta, **kwargs)

    
#pgc = PGCenter()
#pgc = PGCenter((1,2), (3,4))
#pgd = PGData()
fname = '/data/io/IoIO/raw/2020-07-15/HD87696-0016_Na_off.fit'
#pgd = PGData(fname)

#ccd = CCDData.read(fname, unit='adu')

pgd = PGData.read(fname)
cpgd = PGCCDData.read(fname)

rpgd = PGData(pgd.data, meta=pgd.meta)

rcpgd = PGCCDData(pgd.data, meta=pgd.meta)

#rname = '/data/io/IoIO/reduced/Calibration/2020-05-15_V_flat.fits'
#pgccd = PGCCDData.read(rname)
#
#with fits.open(rname) as hdul:
#    hdul_pgccd = PGCCDData(hdul)
#    print(hdul.fileinfo(0))
#    print(len(hdul))

dspgd = PGData.read(fname, desired_center=(100,100))
cdspgd = PGCCDData.read(fname, desired_center=(100,100))

#test = MPGCCDData.read(fname)

print('done')
