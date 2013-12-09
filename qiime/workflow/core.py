#!/usr/bin/env python

class Workflow(object):
    """Arbitrary worflow support structure"""
    def __init__(self, ShortCircuit=True, **kwargs):
        """Build thy self

        ShortCiruit : if True, enables ignoring function groups when a given
            item has failed

        kwargs are stored as self.Options. Support for arbitrary Stats is 
        implicit
        """
        self.Options = kwargs
        self.Stats = defaultdict(int)
        self.ShortCircuit = ShortCircuit
        
    def _construct_iterator(self, **kwargs):
        """Define the central iterator"""
        raise NotImplementedError("Must be implemented")

    def _assign_function_groups(self, **kwargs):
        """Determine what function groups will be used

        A function group is simply a function that subsequently calls the
        methods of interested. For instance, you may have a _process_seqs 
        function group, that then calls _check_length, _split_sequence, etc.
        """
        raise NotImplementedError("Must be implemented")

    def _initialize_item_state(self, item):
        """Initialie the per-item state in self"""
        raise NotImplementedError("Must be implemented")

    def __call__(self, success_callback, failed_callback, **kwargs):
        """Operate on all the data

        success_callback : method to call on a successful item prior to 
            yielding
        failed_callback : method to call on a failed item prior to yielding
        kwargs : these will get passed to the iterator constructor and to the
            the method that determines the function groups
        """
        gen = self._construct_iterator(**kwargs)
        function_groups = self._assign_function_groups(**kwargs)

        for item in gen:
            self._initialize_item_state(item)
            
            for f in function_groups:
                f(item)

            if self.Failed:
                yield failed_callback(self.FinalState)
            else:
                yield success_callback(self.FinalState)

    def requires(self, f, IsValid=True, Option=None, Values=None):
        """Decorator that executes a function if requirements are met

        f : the decorated function
        IsValid : execute the function if self.Failed is False
        Option : a required option
        Values : required values associated with an option
        """
        if not isinstance(values, set):
            if isinstance(values, Iterable)
                values = set(values)
            else:
                values = set([values])
                
        def decorated(self, item):
            if IsValid and self.Failed:
                return
                
            if self.Options[option] in values:
                f(item)
