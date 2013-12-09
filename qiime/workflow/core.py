#!/usr/bin/env python

from collections import Iterable, defaultdict

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
        self.Failed = False
        self.FinalState = None

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

    def __call__(self, success_callback=None, failed_callback=None, **kwargs):
        """Operate on all the data

        success_callback : method to call on a successful item prior to 
            yielding
        failed_callback : method to call on a failed item prior to yielding
        kwargs : these will get passed to the iterator constructor and to the
            the method that determines the function groups
        """
        if success_callback is None:
            success_callback = lambda x: x

        gen = self._construct_iterator(**kwargs)
        function_groups = self._assign_function_groups(**kwargs)

        for item in gen:
            self.Failed = False
            self._initialize_item_state(item)
            
            for f in function_groups:
                f(item)

            if self.Failed and failed_callback is not None:
                yield failed_callback(self.FinalState)
            else:
                yield success_callback(self.FinalState)

class requires(object):
    """Decorator that executes a function if requirements are met"""
    def __init__(self, IsValid=True, Option=None, Values=None):
        """
        f : the decorated function
        IsValid : execute the function if self.Failed is False
        Option : a required option
        Values : required values associated with an option
        """
        # self here is the requires object
        self.IsValid = IsValid
        self.Option = Option
        self.Values = Values

        if not isinstance(self.Values, set):
            if isinstance(self.Values, Iterable):
                self.Values = set(self.Values)
            else:
                self.Values = set([self.Values])
    
    def __call__(outer_self, f):
        # outer_self is the requires object
        # self is expected to be a Workflow object
        def decorated_with_option(self, *args, **kwargs):
            if outer_self.IsValid and (self.Failed and self.ShortCircuit):
                return
            
            opt = self.Options.get(outer_self.Option, 'MISSING OPTION')
            if opt != 'MISSING OPTION' and opt in outer_self.Values:
                f(self, *args, **kwargs)
        
        def decorated_without_option(self, *args, **kwargs):
            if outer_self.IsValid and (self.Failed and self.ShortCircuit):
                return

            f(self, *args, **kwargs)

        if outer_self.Option is not None:
            return decorated_with_option
        else:
            return decorated_without_option
