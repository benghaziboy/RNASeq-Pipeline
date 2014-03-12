class DefaultList(list):

    def __init__(self, *args, **kwargs):
        list.__init__(self, *args)
        self.default = kwargs.get('default', None)

    def __getitem__(self, key):
        # retrieving an item does not expand the list
        if isinstance(key, slice):
            return [self[elt] for elt in range(key.start, key.stop, key.step)]
        else:
            try:
                return list.__getitem__(self, key)
            except IndexError:
                return self.default

    def __setitem__(self, key, value):
        # setting an item may expand the list
        try:
            list.__setitem__(self, key, value)
        except IndexError:
            self.extend([self.default] * (key - len(self)))
            self.append(value)
