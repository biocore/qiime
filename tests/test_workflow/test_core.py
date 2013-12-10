#!/usr/bin/env python

from itertools import izip
from qiime.workflow.core import Workflow, requires
from cogent.util.unit_test import TestCase, main

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
    def _assign_function_groups(self, **kwargs):
        groups = []
        if 'A' in kwargs:
            groups.append(self.groupA)
        if 'B' in kwargs:
            groups.append(self.groupB)
        if 'C' in kwargs:
            groups.append(self.groupC)
        return groups

    def groupA(self, item):
        self.methodA1(item)
        self.methodA2(item)

    def groupB(self, item):
        self.methodB1(item)
        self.methodB2(item)
   
    @requires(IsValid=True)
    def groupC(self, item):
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

    @requires(IsValid=True)
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
        self.obj_short = MockWorkflow(**{'A1':'foo', 'xyz':10,'C2':2})
        self.obj_noshort = MockWorkflow(ShortCircuit=False, **{'A1':'foo', 
                                                             'xyz':10,'C2':2})

    def test_init(self):
        self.assertEqual(self.obj_short.Options, {'A1':'foo', 'xyz':10, 'C2':2})
        self.assertEqual(self.obj_short.Stats, {})
        self.assertTrue(self.obj_short.ShortCircuit)
        self.assertEqual(self.obj_noshort.Options, {'A1':'foo', 'xyz':10, 
                                                   'C2':2})
        self.assertEqual(self.obj_noshort.Stats, {})
        self.assertFalse(self.obj_noshort.ShortCircuit)

    def test_assign_function_groups(self):
        exp_None = []
        exp_AB = [self.obj_short.groupA, self.obj_short.groupB]

        obs_None = self.obj_short._assign_function_groups()
        obs_AB = self.obj_short._assign_function_groups(**{'A':None, 'B':None})
        
        self.assertEqual(obs_None, exp_None)
        self.assertEqual(obs_AB, exp_AB)

    def test_call_AC_no_fail(self):
        single_iter = construct_iterator(**{'iter_x':[1,2,3,4,5]})
        sf = lambda x: x.FinalState # success function
        
        exp_stats = {'A1':5, 'A2':5, 'C1':5, 'C2':5}
        exp_result = [('C2',1), ('C2',2), ('C2',3), ('C2',4), ('C2', 5)]

        kwargs = {'A':None, 'C':None}
        obs_result = list(self.obj_short(single_iter, sf, None, **kwargs))

        self.assertEqual(obs_result, exp_result)
        self.assertEqual(self.obj_short.Stats, exp_stats)

    def test_call_AC_fail(self):
        single_iter = construct_iterator(**{'iter_x':[1,2,'fail A2',4,5]})
        sf = lambda x: x.FinalState # success function
        ff = lambda x: x.FinalState # failed function
        
        exp_stats = {'A1':5, 'A2':5, 'C1':4, 'C2':4}

        kwargs = {'A':None, 'C':None, 'C2':1}

        # pass in a failed callback to capture the result, and pause execution
        gen = self.obj_short(single_iter, sf, ff, **kwargs)

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
        
        exp_stats = {'A1':5, 'A2':5, 'C1':5, 'C2':5}

        kwargs = {'A':None, 'C':None}

        # pass in a failed callback to capture the result, and pause execution
        gen = self.obj_noshort(single_iter, sf, ff, **kwargs)

        r1 = gen.next()
        self.assertEqual(r1, ('C2', 1))
        self.assertFalse(self.obj_noshort.Failed)

        r2 = gen.next()
        self.assertEqual(r2, ('C2', 2))
        self.assertFalse(self.obj_noshort.Failed)
        
        r3 = gen.next()
        self.assertEqual(self.obj_noshort.FinalState, ('C2', 'fail A2'))
        self.assertTrue(self.obj_noshort.Failed)

        r4 = gen.next()
        self.assertEqual(r4, ('C2', 4))
        self.assertFalse(self.obj_noshort.Failed)
        
        r5 = gen.next()
        self.assertEqual(r5, ('C2', 5))
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
        
if __name__ == '__main__':
    main()
