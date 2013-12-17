#!/usr/bin/env python

"""Perform multiple method calls, determined at runtime, on independent items

Construct arbitrarily complex workflows in which the specific methods run are
determined at runtime. These methods are applied to items that are assumed to
be independent.

As an example:

class MyWorkflow(Workflow):
    @priority(100)
    @no_requirements
    def wf_mul(self, item):
        self.FinalState *= item

    @priority(10)
    @requires(Option='add_value')
    def wf_add(self, item):
        self.FinalState += item
    
    @requires(Option='sub_value', Values=[1,5,10])
    def wf_sub(self, item):
        self.FinalState -= item
        self.FinalState -= self.Options['sub_value']

    @priority(1000)
    @requires(IsValid=False)
    def wf_init(self, item):
        self.FinalState = item

# (i * i) + i - i - 5
wf = MyWorkflow(Options={'add_value':None, 'sub_value':5})
gen = (i for i in range(10))
for i in wf(gen):
    print i

# (i * i) - i - 10
wf = MyWorkflow(Options={'sub_value':10})
gen = (i for i in range(10))
for i in wf(gen):
    print i

# (i * i)
wf = MyWorkflow()
gen = (i for i in range(10))
for i in wf(gen):
    print i
"""

from itertools import chain
from functools import update_wrapper
from collections import Iterable, defaultdict

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2013, The QIIME Project"
__credits__ = ["Daniel McDonald", "Tony Walters"]
__license__ = "BSD" # NOTE, this script does _not_ import GPL code
__version__ = "1.7.0-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"

# thank you Flask project...
_missing = object()
_executed = object()

class Exists(object):
    def __contains__(self, item):
        return True
option_exists = Exists()

class Workflow(object):
    """Arbitrary worflow support structure"""

    def __init__(self, ShortCircuit=True, Debug=True, Options=None, **kwargs):
        """Build thy self

        ShortCiruit : if True, enables ignoring function groups when a given
            item has failed
        Debug : Enable debug mode
        Options : runtime options, {'option':values}
        kwargs : Additional arguments will be added to self

        All workflow methods (i.e., those starting with "wk_") must be decorated
        by either "no_requirements" or "requires". This ensures that the methods
        support the automatic workflow determination mechanism.
        """
        if Options is None:
            self.Options = {}
        else:
            self.Options = Options

        self.Stats = defaultdict(int)
        self.ShortCircuit = ShortCircuit
        self.Failed = False
        self.FinalState = None
        
        for k,v in kwargs.iteritems():
            if hasattr(self, k):
                raise AttributeError("%s exists in self!" % k)
            setattr(self, k, v)

        for f in self._all_wf_methods():
            if not hasattr(f, '__workflowtag__'):
                raise AttributeError("%s isn't a workflow method!" % f.__name__)

        self._stage_state()
        self._sanity_check()

    def _stage_state(self):
        """Stage any additional data necessary for the workflow
        
        This does not need to be overloaded
        """
        pass

    def _sanity_check(self):
        """Perform a sanity check on self"""
        raise NotImplementedError("Must implement a sanity check!")

    def _all_wf_methods(self, default_priority=0):
        """Get all workflow methods
        
        Methods are sorted by priority
        """
        methods = [getattr(self, f) for f in dir(self) if f.startswith('wf_')]
        key = lambda x: getattr(x, 'Priority', default_priority)
        return sorted(methods, key=key, reverse=True)

    def _get_workflow(self, it):
        """Get the methods executed, sorted by priority"""
        # save state
        shortcircuit_state = self.ShortCircuit
        self.ShortCircuit = False
        stats = self.Stats.copy()

        peek = it.next()
        executed = [f for f in self._all_wf_methods() if f(peek) is _executed]

        # restore state
        self.ShortCircuit = shortcircuit_state
        self.Stats = stats
        generator_reset = chain([peek], it)

        return generator_reset, executed

    def __call__(self, it, success_callback=None, fail_callback=None):
        """Operate on all the data

        it : an iterator
        success_callback : method to call on a successful item prior to 
            yielding
        fail_callback : method to call on a failed item prior to yielding
        """
        if success_callback is None:
            success_callback = lambda x: x.FinalState

        it, workflow = self._get_workflow(it)
        
        for item in it:
            self.Failed = False

            for f in workflow:
                f(item)

            if self.Failed and fail_callback is not None:
                yield fail_callback(self)
            else:
                yield success_callback(self)
            
    @staticmethod
    def tagFunction(f):
        setattr(f, '__workflowtag__', None)

def no_requirements(f):
    def decorated(self, *args, **kwargs):
        f(self, *args, **kwargs)
        return _executed
    Workflow.tagFunction(decorated)
    return update_wrapper(decorated, f)

class requires(object):
    """Decorator that executes a function if requirements are met"""
    def __init__(self, IsValid=True, Option=None, Values=_missing):
        """
        IsValid : execute the function if self.Failed is False
        Option : a required option
        Values : required values associated with an option
        """
        # self here is the requires object
        self.IsValid = IsValid
        self.Option = Option

        if Values is _missing:
            self.Values = option_exists
        elif not isinstance(Values, set):
            if isinstance(Values, Iterable):
                self.Values = set(Values)
            else:
                self.Values = set([Values])
        else:
            self.Values = Values

    def doShortCircuit(self, wrapped):
        if self.IsValid and (wrapped.Failed and wrapped.ShortCircuit):
            return True
        else:
            return False

    def __call__(self, f):
        """Wrap a function

        f : the function to wrap
        """
        def decorated_with_option(dec_self, *args, **kwargs):
            """A decorated function that has an option to validate

            dec_self : this is "self" for the decorated function
            """
            if self.doShortCircuit(dec_self):
                return

            s_opt = self.Option
            ds_opts = dec_self.Options

            if s_opt in ds_opts and ds_opts[s_opt] in self.Values:
                f(dec_self, *args, **kwargs)
                return _executed

        def decorated_without_option(dec_self, *args, **kwargs):
            """A decorated function that does not have an option to validate

            dec_self : this is "self" for the decorated function
            """
            if self.doShortCircuit(dec_self):    
                return

            f(dec_self, *args, **kwargs)
            return _executed

        Workflow.tagFunction(decorated_with_option)
        Workflow.tagFunction(decorated_without_option)

        if self.Option is None:
            return update_wrapper(decorated_without_option, f)
        else:
            return update_wrapper(decorated_with_option, f)

class priority(object):
    """Sets a function priority"""
    def __init__(self, Priority):
        self.Priority = Priority

    def __call__(self, f):
        f.Priority = self.Priority
        return f
