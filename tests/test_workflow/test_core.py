#!/usr/bin/env python

from itertools import izip
from qiime.workflow.core import (Workflow, requires, priority,
        no_requirements)
from unittest import TestCase, main

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2013, The QIIME Project"
__credits__ = ["Daniel McDonald"]
__license__ = "BSD" # NOTE, does not import any GPL code
__version__ = "1.7.0-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"

def construct_iterator(**kwargs):
    """make an iterator for testing purposes"""
    to_gen = []
    for k in sorted(kwargs):
        if k.startswith('iter'):
            to_gen.append(kwargs[k])
    if len(to_gen) == 1:
        return (x for x in to_gen[0])
    else:
        return izip(*to_gen)

class MockWorkflow(Workflow):
    @priority(90)
    @requires(Option='A', Values=True)
    def wf_groupA(self, item):
        self.methodA1(item)
        self.methodA2(item)

    @requires(Option='B', Values=True)
    def wf_groupB(self, item):
        self.methodB1(item)
        self.methodB2(item)
   
    @priority(10)
    @requires(Option='C', Values=True)
    def wf_groupC(self, item):
        self.methodC1(item)
        self.methodC2(item)

    @requires(IsValid=False) # always execute
    def methodA1(self, item):
        name = 'A1'
        self.Stats[name] += 1
        if item == 'fail %s' % name:
            self.Failed = True
        self.FinalState = (name, item)

    def methodA2(self, item):
        name = 'A2'
        self.Stats[name] += 1
        if item == 'fail %s' % name:
            self.Failed = True
        self.FinalState = (name, item)

    @requires(IsValid=False)
    def methodB1(self, item):
        name = 'B1'
        self.Stats[name] += 1
        if item == 'fail %s' % name:
            self.Failed = True
            self.FinalState = 'failed'
        else:
            self.FinalState = (name, item)

    @requires(Option='foo', Values=[1,2,3])
    def methodB2(self, item):
        name = 'B2'
        self.Stats[name] += 1
        if item == 'fail %s' % name:
            self.Failed = True
            self.FinalState = 'failed'
        else:
            self.FinalState = (name, item)

    @no_requirements
    def methodC1(self, item):
        name = 'C1'
        self.Stats[name] += 1
        if item == 'fail %s' % name:
            self.Failed = True
        self.FinalState = (name, item)

    @requires(IsValid=True, Option='C2', Values=[1,2,3])
    def methodC2(self, item):
        name = 'C2'
        self.Stats[name] += 1
        if item == 'fail %s' % name:
            self.Failed = True
        self.FinalState = (name, item)

class WorkflowTests(TestCase):
    def setUp(self):
        self.obj_short = MockWorkflow(**{'A':True, 'C':True})
        self.obj_noshort = MockWorkflow(ShortCircuit=False, **{'A':True, 
                                                               'C':True})

    def test_untagged_wf_method(self):
        class WFTest(Workflow):
            @no_requirements
            def wf_1(self):
                pass
            def wf_2(self):
                pass

        with self.assertRaises(AttributeError):
            _ = WFTest()

    def test_get_workflow(self):
        gen = single_iter = construct_iterator(**{'iter_x':[1,2,3,4,5]})
        exp_wf = [self.obj_short.wf_groupA, self.obj_short.wf_groupC]
        obs_gen, obs_wf = self.obj_short._get_workflow(gen)

        self.assertEqual(obs_wf, exp_wf)
        self.assertEqual(list(obs_gen), [1,2,3,4,5])
        
        self.assertEqual(self.obj_short.Stats, {})
        self.assertTrue(self.obj_short.ShortCircuit)
    
    def test_init(self):
        self.assertEqual(self.obj_short.Options, {'A':True, 'C':True})
        self.assertEqual(self.obj_short.Stats, {})
        self.assertTrue(self.obj_short.ShortCircuit)
        self.assertEqual(self.obj_noshort.Options, {'A':True, 'C':True})
        self.assertEqual(self.obj_noshort.Stats, {})
        self.assertFalse(self.obj_noshort.ShortCircuit)

    def test_all_wf_methods(self):
        # note on priority: groupA:90, groupC:10, groupB:0 (default)
        exp = [self.obj_short.wf_groupA, self.obj_short.wf_groupC,
               self.obj_short.wf_groupB]
        obs = self.obj_short._all_wf_methods()
        self.assertEqual(obs, exp)

    def test_call_AC_no_fail(self):
        single_iter = construct_iterator(**{'iter_x':[1,2,3,4,5]})
        sf = lambda x: x.FinalState # success function
        
        exp_stats = {'A1':5, 'A2':5, 'C1':5}
        # C2 isn't executed as its requirements aren't met in the Options
        exp_result = [('C1',1), ('C1',2), ('C1',3), ('C1',4), ('C1', 5)]

        obs_result = list(self.obj_short(single_iter, sf, None))

        self.assertEqual(obs_result, exp_result)
        self.assertEqual(self.obj_short.Stats, exp_stats)

    def test_call_AC_fail(self):
        single_iter = construct_iterator(**{'iter_x':[1,2,'fail A2',4,5]})
        sf = lambda x: x.FinalState # success function
        ff = lambda x: x.FinalState # failed function
        
        exp_stats = {'A1':5, 'A2':5, 'C1':4, 'C2':4}

        self.obj_short.Options['C2'] = 1
        # pass in a failed callback to capture the result, and pause execution
        gen = self.obj_short(single_iter, sf, ff)

        r1 = gen.next()
        self.assertEqual(r1, ('C2', 1))
        self.assertFalse(self.obj_short.Failed)

        r2 = gen.next()
        self.assertEqual(r2, ('C2', 2))
        self.assertFalse(self.obj_short.Failed)
        
        r3 = gen.next()
        self.assertEqual(self.obj_short.FinalState, ('A2', 'fail A2'))
        self.assertTrue(self.obj_short.Failed)
        self.assertEqual(r3, ('A2', 'fail A2'))

        r4 = gen.next()
        self.assertEqual(r4, ('C2', 4))
        self.assertFalse(self.obj_short.Failed)
        
        r5 = gen.next()
        self.assertEqual(r5, ('C2', 5))
        self.assertFalse(self.obj_short.Failed)

        self.assertEqual(self.obj_short.Stats, exp_stats)

    def test_call_AC_fail_noshort(self):
        single_iter = construct_iterator(**{'iter_x':[1,2,'fail A2',4,5]})
        sf = lambda x: x.FinalState # success function
        ff = lambda x: x.FinalState # failed function
        
        exp_stats = {'A1':5, 'A2':5, 'C1':5}

        # pass in a failed callback to capture the result, and pause execution
        gen = self.obj_noshort(single_iter, sf, ff)

        r1 = gen.next()
        self.assertEqual(r1, ('C1', 1))
        self.assertFalse(self.obj_noshort.Failed)

        r2 = gen.next()
        self.assertEqual(r2, ('C1', 2))
        self.assertFalse(self.obj_noshort.Failed)
        
        r3 = gen.next()
        self.assertEqual(self.obj_noshort.FinalState, ('C1', 'fail A2'))
        self.assertTrue(self.obj_noshort.Failed)

        r4 = gen.next()
        self.assertEqual(r4, ('C1', 4))
        self.assertFalse(self.obj_noshort.Failed)
        
        r5 = gen.next()
        self.assertEqual(r5, ('C1', 5))
        self.assertFalse(self.obj_noshort.Failed)

        self.assertEqual(self.obj_noshort.Stats, exp_stats)

class RequiresTests(TestCase):
    def test_methodb1(self):
        obj = MockWorkflow()
        obj.methodB1('test')
        self.assertEqual(obj.FinalState, ('B1', 'test'))
        self.assertFalse(obj.Failed)
        
        # methodb1 executes regardless of if self.Failed
        obj.Failed = True
        obj.methodB1('test 2')
        self.assertEqual(obj.FinalState, ('B1', 'test 2'))

        obj.Failed = False
        obj.methodB1('fail B1')
        self.assertEqual(obj.FinalState, 'failed')
       
        self.assertEqual(obj.Stats, {'B1':3})

    def test_methodb2_accept(self):
        # methodb2 is setup to be valid when foo is in [1,2,3], make sure we
        # can execute
        obj = MockWorkflow(**{'foo':1})
        obj.methodB2('test')
        self.assertEqual(obj.FinalState, ('B2', 'test'))
        self.assertEqual(obj.Stats, {'B2':1})

        # methodb2 will not execute if self.Failed
        obj.Failed = True
        obj.methodB2('test 2')
        self.assertEqual(obj.FinalState, ('B2', 'test'))
        self.assertEqual(obj.Stats, {'B2':1})

    def test_methodb2_ignore(self):
        # methodb2 is setup to be valid when foo is in [1, 2, 3], make sure
        # we do not execute
        obj = MockWorkflow(**{'foo':'bar'})
        obj.methodB2('test')
        self.assertEqual(obj.FinalState, None)
        self.assertEqual(obj.Stats, {})

class PriorityTests(TestCase):
    def test_dec(self):
        @priority(10)
        def foo(x,y,z):
            """doc check"""
            return x+y+z

        self.assertEqual(foo.Priority, 10)
        self.assertEqual(foo.__name__, 'foo')
        self.assertEqual(foo.__doc__, 'doc check')

if __name__ == '__main__':
    main()
