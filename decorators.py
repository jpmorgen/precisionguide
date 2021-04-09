"""Decorators for precisionguide system"""

import numpy as np

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
        . tuple: shape of resulting np.array must match shape_check or
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

