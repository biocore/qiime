#!/usr/bin/env python
"""copy of cogent.core.tree to provide newer getSubTree() and parse_newick
"""
from numpy import zeros, argsort
from copy import deepcopy
import re
from cogent.util.transform import comb
from cogent.maths.stats.test import correlation
from operator import or_
from cogent.util.misc import InverseDict
from random import shuffle

__author__ = "Gavin Huttley, Peter Maxwell and Rob Knight"
__copyright__ = "Copyright 2007-2011, The Cogent Project"
__credits__ = ["Gavin Huttley", "Peter Maxwell", "Rob Knight",
                    "Andrew Butterfield", "Catherine Lozupone", "Micah Hamady",
                    "Jeremy Widmann", "Zongzhi Liu", "Daniel McDonald",
                    "Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.6.0.dev"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

def distance_from_r_squared(m1, m2):
    """Estimates distance as 1-r^2: no correl = max distance"""
    return 1 - (correlation(m1.flat, m2.flat)[0])**2

def distance_from_r(m1, m2):
    """Estimates distance as (1-r)/2: neg correl = max distance"""
    return (1-correlation(m1.flat, m2.flat)[0])/2

class TreeError(Exception):
    pass

class TreeNode(object):
    """Store information about a tree node. Mutable.
    
    Parameters:
        Name: label for the node, assumed to be unique.
        Children: list of the node's children.
        Params: dict containing arbitrary parameters for the node.
        NameLoaded: ?
    """
    _exclude_from_copy = dict.fromkeys(['_parent','Children'])
    
    def __init__(self, Name=None, Children=None, Parent=None, Params=None, \
            NameLoaded=True, **kwargs):
        """Returns new TreeNode object."""
        self.Name = Name
        self.NameLoaded = NameLoaded
        if Params is None:
            Params = {}
        self.params = Params
        self.Children = []
        if Children is not None:
            self.extend(Children)
        self._parent = Parent
        if (Parent is not None) and not (self in Parent.Children):
            Parent.append(self)

### built-in methods and list interface support
    def __repr__(self):
        """Returns reconstructable string representation of tree.
        
        WARNING: Does not currently set the class to the right type.
        """
        return 'Tree("%s")' % self.getNewick()
    
    def __str__(self):
        """Returns Newick-format string representation of tree."""
        return self.getNewick()
    
    def compareName(self, other):
        """Compares TreeNode by name"""
        if self is other:
            return 0
        try:
            return cmp(self.Name, other.Name)
        except AttributeError:
            return cmp(type(self), type(other))

    def compareByNames(self, other):
        """Equality test for trees by name"""
        # if they are the same object then they must be the same tree...
        if self is other:
            return True
        self_names = self.getNodeNames()
        other_names = other.getNodeNames()
        self_names.sort()
        other_names.sort()
        return self_names == other_names

    def _to_self_child(self, i):
        """Converts i to self's type, with self as its parent.
        
        Cleans up refs from i's original parent, but doesn't give self ref to i.
        """
        c = self.__class__
        if isinstance(i, c):
            if i._parent not in (None, self):
                i._parent.Children.remove(i)
        else:
            i = c(i)
        i._parent = self
        return i
    
    def append(self, i):
        """Appends i to self.Children, in-place, cleaning up refs."""
        self.Children.append(self._to_self_child(i))
    
    def extend(self, items):
        """Extends self.Children by items, in-place, cleaning up refs."""
        self.Children.extend(map(self._to_self_child, items))
    
    def insert(self, index, i):
        """Inserts an item at specified position in self.Children."""
        self.Children.insert(index, self._to_self_child(i))
    
    def pop(self, index=-1):
        """Returns and deletes child of self at index (default: -1)"""
        result = self.Children.pop(index)
        result._parent = None
        return result
    
    def remove(self, target):
        """Removes node by name instead of identity.
        
        Returns True if node was present, False otherwise.
        """
        if isinstance(target, TreeNode):
            target = target.Name
        for (i, curr_node) in enumerate(self.Children):
            if curr_node.Name == target:
                self.removeNode(curr_node)
                return True
        return False
    
    def __getitem__(self, i):
        """Node delegates slicing to Children; faster to access them
        directly."""
        return self.Children[i]
    
    def __setitem__(self, i, val):
        """Node[i] = x sets the corresponding item in Children."""
        curr = self.Children[i]
        if isinstance(i, slice):
            for c in curr:
                c._parent = None
            coerced_val = map(self._to_self_child, val)
            self.Children[i] = coerced_val[:]
        else:   #assume we got a single index
            curr._parent = None
            coerced_val = self._to_self_child(val)
            self.Children[i] = coerced_val
    
    def __delitem__(self, i):
        """del node[i] deletes index or slice from self.Children."""
        curr = self.Children[i]
        if isinstance(i, slice):
            for c in curr:
                c._parent = None
        else:
            curr._parent = None
        del self.Children[i]
    
    def __iter__(self):
        """Node iter iterates over the Children."""
        return iter(self.Children)
    
    def __len__(self):
        """Node len returns number of children."""
        return len(self.Children)
    
    #support for copy module
    def copyRecursive(self, memo=None, _nil=[], constructor='ignored'):
        """Returns copy of self's structure, including shallow copy of attrs.
        
        constructor is ignored; required to support old tree unit tests.
        """
        result = self.__class__()
        efc = self._exclude_from_copy
        for k, v in self.__dict__.items():
            if k not in efc:  #avoid infinite recursion
                result.__dict__[k] = deepcopy(self.__dict__[k])
        for c in self:
            result.append(c.copy())
        return result
    
    def copy(self, memo=None, _nil=[], constructor='ignored'):
        """Returns a copy of self using an iterative approach"""
        def __copy_node(n):
            result = n.__class__()
            efc = n._exclude_from_copy
            for k,v in n.__dict__.items():
                if k not in efc:
                    result.__dict__[k] = deepcopy(n.__dict__[k])
            return result

        root = __copy_node(self)
        nodes_stack = [[root, self, len(self.Children)]]

        while nodes_stack:
            #check the top node, any children left unvisited?
            top = nodes_stack[-1]
            new_top_node, old_top_node, unvisited_children = top

            if unvisited_children:
                top[2] -= 1
                old_child = old_top_node.Children[-unvisited_children]
                new_child = __copy_node(old_child)
                new_top_node.append(new_child)
                nodes_stack.append([new_child, old_child, \
                                    len(old_child.Children)])
            else:  #no unvisited children
                nodes_stack.pop()
        return root

    __deepcopy__ = deepcopy = copy
   
    def copyTopology(self, constructor=None):
        """Copies only the topology and labels of a tree, not any extra data.
        
        Useful when you want another copy of the tree with the same structure
        and labels, but want to e.g. assign different branch lengths and
        environments. Does not use deepcopy from the copy module, so _much_
        faster than the copy() method.
        """
        if constructor is None:
            constructor = self.__class__
        children = [c.copyTopology(constructor) for c in self.Children]
        return constructor(Name=self.Name[:], Children=children)
    
    #support for basic tree operations -- finding objects and moving in the tree
    def _get_parent(self):
        """Accessor for parent.
        
        If using an algorithm that accesses Parent a lot, it will be much
        faster to access self._parent directly, but don't do it if mutating
        self._parent! (or, if you must, remember to clean up the refs).
        """
        return self._parent
    
    def _set_parent(self, Parent):
        """Mutator for parent: cleans up refs in old parent."""
        if self._parent is not None:
            self._parent.removeNode(self)
        self._parent = Parent
        if (Parent is not None) and (not self in Parent.Children):
            Parent.Children.append(self)
    
    Parent = property(_get_parent, _set_parent)
    
    def indexInParent(self):
        """Returns index of self in parent."""
        return self._parent.Children.index(self)
    
    def isTip(self):
        """Returns True if the current node is a tip, i.e. has no children."""
        return not self.Children
    
    def isRoot(self):
        """Returns True if the current is a root, i.e. has no parent."""
        return self._parent is None
   
    def traverse(self, self_before=True, self_after=False, include_self=True):
        """Returns iterator over descendants. Iterative: safe for large trees.

        self_before includes each node before its descendants if True.
        self_after includes each node after its descendants if True.
        include_self includes the initial node if True.
                    
        self_before and self_after are independent. If neither is True, only
        terminal nodes will be returned.
                    
        Note that if self is terminal, it will only be included once even if
        self_before and self_after are both True.
            
        This is a depth-first traversal. Since the trees are not binary,
        preorder and postorder traversals are possible, but inorder traversals
        would depend on the data in the tree and are not handled here.
        """
        if self_before:
            if self_after:
                return self.pre_and_postorder(include_self=include_self)
            else:
                return self.preorder(include_self=include_self)
        else:
            if self_after:
                return self.postorder(include_self=include_self)
            else:
                return self.tips(include_self=include_self)

    def levelorder(self, include_self=True):
        """Performs levelorder iteration over tree"""
        queue = [self]
        while queue:
            curr = queue.pop(0)
            if include_self or (curr is not self):
                yield curr
            if curr.Children:
                queue.extend(curr.Children)

    def preorder(self, include_self=True):
        """Performs preorder iteration over tree."""
        stack = [self]
        while stack:
            curr = stack.pop()
            if include_self or (curr is not self):
                yield curr
            if curr.Children:
                stack.extend(curr.Children[::-1])   #20% faster than reversed    
    def postorder(self, include_self=True):
        """Performs postorder iteration over tree.
        
        This is somewhat inelegant compared to saving the node and its index
        on the stack, but is 30% faster in the average case and 3x faster in
        the worst case (for a comb tree).

        Zongzhi Liu's slower but more compact version is:

        def postorder_zongzhi(self):
            stack = [[self, 0]]
            while stack:
                curr, child_idx = stack[-1]
                if child_idx < len(curr.Children):
                    stack[-1][1] += 1
                    stack.append([curr.Children[child_idx], 0])
                else:
                    yield stack.pop()[0]
        """
        child_index_stack = [0]
        curr = self
        curr_children = self.Children
        curr_children_len = len(curr_children)
        while 1:
            curr_index = child_index_stack[-1]
            #if there are children left, process them
            if curr_index < curr_children_len:
                curr_child = curr_children[curr_index]
                #if the current child has children, go there
                if curr_child.Children:
                    child_index_stack.append(0)
                    curr = curr_child
                    curr_children = curr.Children
                    curr_children_len = len(curr_children)
                    curr_index = 0
                #otherwise, yield that child
                else:
                    yield curr_child
                    child_index_stack[-1] += 1
            #if there are no children left, return self, and move to
            #self's parent
            else:
                if include_self or (curr is not self):
                    yield curr
                if curr is self:
                    break
                curr = curr.Parent
                curr_children = curr.Children
                curr_children_len = len(curr_children)
                child_index_stack.pop()
                child_index_stack[-1] += 1

    def pre_and_postorder(self, include_self=True):
        """Performs iteration over tree, visiting node before and after."""
        #handle simple case first
        if not self.Children:
            if include_self:
                yield self
            raise StopIteration
        child_index_stack = [0]
        curr = self
        curr_children = self.Children
        while 1:
            curr_index = child_index_stack[-1]
            if not curr_index:
                if include_self or (curr is not self):
                    yield curr
            #if there are children left, process them
            if curr_index < len(curr_children):
                curr_child = curr_children[curr_index]
                #if the current child has children, go there
                if curr_child.Children:
                    child_index_stack.append(0)
                    curr = curr_child
                    curr_children = curr.Children
                    curr_index = 0
                #otherwise, yield that child
                else:
                    yield curr_child
                    child_index_stack[-1] += 1
            #if there are no children left, return self, and move to
            #self's parent
            else:
                if include_self or (curr is not self):
                    yield curr
                if curr is self:
                    break
                curr = curr.Parent
                curr_children = curr.Children
                child_index_stack.pop()
                child_index_stack[-1] += 1
 
    def traverse_recursive(self, self_before=True, self_after=False, \
        include_self=True):
        """Returns iterator over descendants. IMPORTANT: read notes below.

        traverse_recursive is slower than traverse, and can lead to stack
        errors. However, you _must_ use traverse_recursive if you plan to
        modify the tree topology as you walk over it (e.g. in post-order),
        because the iterative methods use their own stack that is not updated
        if you alter the tree.
        
        self_before includes each node before its descendants if True.
        self_after includes each node after its descendants if True.
        include_self includes the initial node if True.
        
        self_before and self_after are independent. If neither is True, only
        terminal nodes will be returned.
        
        Note that if self is terminal, it will only be included once even if
        self_before and self_after are both True.
        
        This is a depth-first traversal. Since the trees are not binary,
        preorder and postorder traversals are possible, but inorder traversals
        would depend on the data in the tree and are not handled here.
        """
        if self.Children:
            if self_before and include_self:
                yield self
            for child in self.Children:
                for i in child.traverse_recursive(self_before, self_after):
                    yield i
            if self_after and include_self:
                yield self
        elif include_self:
            yield self
    
    def ancestors(self):
        """Returns all ancestors back to the root. Dynamically calculated."""
        result = []
        curr = self._parent
        while curr is not None:
            result.append(curr)
            curr = curr._parent
        return result
    
    def root(self):
        """Returns root of the tree self is in. Dynamically calculated."""
        curr = self
        while curr._parent is not None:
            curr = curr._parent
        return curr
    
    def isroot(self):
        """Returns True if root of a tree, i.e. no parent."""
        return self._parent is None
    
    def siblings(self):
        """Returns all nodes that are children of the same parent as self.
        
        Note: excludes self from the list. Dynamically calculated.
        """
        if self._parent is None:
            return []
        result = self._parent.Children[:]
        result.remove(self)
        return result
    
    def iterTips(self, include_self=False):
        """Iterates over tips descended from self, [] if self is a tip."""
        #bail out in easy case
        if not self.Children:
            if include_self:
                yield self
            raise StopIteration
        #use stack-based method: robust to large trees
        stack = [self]
        while stack:
            curr = stack.pop()
            if curr.Children:
                stack.extend(curr.Children[::-1])   #20% faster than reversed
            else:
                yield curr 

    def tips(self, include_self=False):
        """Returns tips descended from self, [] if self is a tip."""
        return list(self.iterTips(include_self=include_self))

    def iterNontips(self, include_self=False):
        """Iterates over nontips descended from self, [] if none.

        include_self, if True (default is False), will return the current
        node as part of the list of nontips if it is a nontip."""
        for n in self.traverse(True, False, include_self):
            if n.Children:
                yield n

    def nontips(self, include_self=False):
        """Returns nontips descended from self."""
        return list(self.iterNontips(include_self=include_self))
    
    def istip(self):
        """Returns True if is tip, i.e. no children."""
        return not self.Children
    
    def tipChildren(self):
        """Returns direct children of self that are tips."""
        return [i for i in self.Children if not i.Children]
    
    def nonTipChildren(self):
        """Returns direct children in self that have descendants."""
        return [i for i in self.Children if i.Children]
    
    def childGroups(self):
        """Returns list containing lists of children sharing a state.
        
        In other words, returns runs of tip and nontip children.
        """
        #bail out in trivial cases of 0 or 1 item
        if not self.Children:
            return []
        if len(self.Children) == 1:
            return [self.Children[0]]
        #otherwise, have to do it properly...
        result = []
        curr = []
        state = None
        for i in self.Children:
            curr_state = bool(i.Children)
            if curr_state == state:
                curr.append(i)
            else:
                if curr:
                    result.append(curr)
                    curr = []
                curr.append(i)
                state = curr_state
        #handle last group
        result.append(curr)
        return result
    
    def lastCommonAncestor(self, other):
        """Finds last common ancestor of self and other, or None.
        
        Always tests by identity.
        """
        my_lineage = set([id(node) for node in [self] + self.ancestors()])
        curr = other
        while curr is not None:
            if id(curr) in my_lineage:
                return curr
            curr = curr._parent
        return None
   
    def lowestCommonAncestor(self, tipnames):
        """Lowest common ancestor for a list of tipnames

        This should be around O(H sqrt(n)), where H is height and n is the
        number of tips passed in.
        """
        if len(tipnames) == 1:
            return self.getNodeMatchingName(tipnames[0])

        tipnames = set(tipnames)
        tips = [tip for tip in self.tips() if tip.Name in tipnames]

        if len(tips) == 0:
            return None

        for t in tips:
            prev = t
            curr = t.Parent

            while curr and not hasattr(curr,'black'):
                setattr(curr,'black',[prev])
                prev = curr
                curr = curr.Parent

            # increase black count, multiple children lead to here
            if curr:
                curr.black.append(prev)

        curr = self
        while len(curr.black) == 1:
            curr = curr.black[0]

        return curr

    lca = lastCommonAncestor #for convenience
    
    #support for more advanced tree operations
    
    def separation(self, other):
        """Returns number of edges separating self and other."""
        #detect trivial case
        if self is other:
            return 0
        #otherwise, check the list of ancestors
        my_ancestors = dict.fromkeys(map(id, [self] + self.ancestors()))
        count = 0
        while other is not None:
            if id(other) in my_ancestors:
                #need to figure out how many steps there were back from self
                curr = self
                while not(curr is None or curr is other):
                    count += 1
                    curr = curr._parent
                return count
            else:
                count += 1
                other = other._parent
        return None
    
    def descendantArray(self, tip_list=None):
        """Returns numpy array with nodes in rows and descendants in columns.
        
        A value of 1 indicates that the decendant is a descendant of that node/
        A value of 0 indicates that it is not
        
        Also returns a list of nodes in the same order as they are listed
        in the array.
        
        tip_list is a list of the names of the tips that will be considered,
        in the order they will appear as columns in the final array. Internal
        nodes will appear as rows in preorder traversal order.
        """
        
        #get a list of internal nodes
        node_list = [node for node in self.traverse() if node.Children]
        node_list.sort()
        
        #get a list of tip names if one is not supplied
        if not tip_list:
            tip_list = [n.Name for n in self.tips()]
            tip_list.sort()
        #make a blank array of the right dimensions to alter
        result = zeros([len(node_list), len(tip_list)])
        #put 1 in the column for each child of each node
        for (i, node) in enumerate(node_list):
            children = [n.Name for n in node.tips()]
            for (j, dec) in enumerate(tip_list):
                if dec in children:
                    result[i,j] = 1
        return result, node_list
    
    def _default_tree_constructor(self):
        return TreeBuilder(constructor=self.__class__).edgeFromEdge
    
    def nameUnnamedNodes(self):
        """sets the Data property of unnamed nodes to an arbitrary value
        
        Internal nodes are often unnamed and so this function assigns a
        value for referencing."""
        #make a list of the names that are already in the tree
        names_in_use = []
        for node in self.traverse():
            if node.Name:
                names_in_use.append(node.Name)
        #assign unique names to the Data property of nodes where Data = None
        name_index = 1
        for node in self.traverse():
            if not node.Name:
                new_name = 'node' + str(name_index)
                #choose a new name if name is already in tree
                while new_name in names_in_use:
                    name_index += 1
                    new_name = 'node' + str(name_index)
                node.Name = new_name
                names_in_use.append(new_name)
                name_index += 1
    
    def makeTreeArray(self, dec_list=None):
        """Makes an array with nodes in rows and descendants in columns.
        
        A value of 1 indicates that the decendant is a descendant of that node/
        A value of 0 indicates that it is not
        
        also returns a list of nodes in the same order as they are listed
        in the array"""
        #get a list of internal nodes
        node_list = [node for node in self.traverse() if node.Children]
        node_list.sort()
        
        #get a list of tips() Name if one is not supplied
        if not dec_list:
            dec_list = [dec.Name for dec in self.tips()]
            dec_list.sort()
        #make a blank array of the right dimensions to alter
        result = zeros((len(node_list), len(dec_list)))
        #put 1 in the column for each child of each node
        for i, node in enumerate(node_list):
            children = [dec.Name for dec in node.tips()]
            for j, dec in enumerate(dec_list):
                if dec in children:
                    result[i,j] = 1
        return result, node_list
    
    def removeDeleted(self,is_deleted):
        """Removes all nodes where is_deleted tests true.
        
        Internal nodes that have no children as a result of removing deleted
        are also removed.
        """
        #Traverse tree
        for node in list(self.traverse(self_before=False,self_after=True)):
            #if node is deleted
            if is_deleted(node):
                #Store current parent
                curr_parent=node.Parent
                #Set current node's parent to None (this deletes node)
                node.Parent=None
                #While there are no chilren at node and not at root
                while (curr_parent is not None) and (not curr_parent.Children):
                    #Save old parent
                    old_parent=curr_parent
                    #Get new parent
                    curr_parent=curr_parent.Parent
                    #remove old node from tree
                    old_parent.Parent=None
    
    def prune(self):
        """Reconstructs correct topology after nodes have been removed.
        
        Internal nodes with only one child will be removed and new connections
        will be made to reflect change.
        """
        #traverse tree to decide nodes to be removed.
        nodes_to_remove = []
        for node in self.traverse():
            if (node.Parent is not None) and (len(node.Children)==1):
                nodes_to_remove.append(node)
        for node in nodes_to_remove:
            #save current parent
            curr_parent=node.Parent
            #save child
            child=node.Children[0]
            #remove current node by setting parent to None
            node.Parent=None
            #Connect child to current node's parent
            child.Parent=curr_parent
    
    def sameShape(self, other):
        """Ignores lengths and order, so trees should be sorted first"""
        if len(self.Children) != len(other.Children):
            return False
        if self.Children:
            for (self_child, other_child) in zip(self.Children, other.Children):
                if not self_child.sameShape(other_child):
                    return False
            return True
        else:
            return self.Name == other.Name
    
    def getNewickRecursive(self, with_distances=False, semicolon=True, \
            escape_name=True):
        """Return the newick string for this edge.
        
        Arguments:
            - with_distances: whether branch lengths are included.
            - semicolon: end tree string with a semicolon
            - escape_name: if any of these characters []'"(),:;_ exist in a
                nodes name, wrap the name in single quotes
        """
        newick = []
        
        subtrees = [child.getNewick(with_distances, semicolon=False)
                for child in self.Children]
        if subtrees:
            newick.append("(%s)" % ",".join(subtrees))
        
        if self.NameLoaded:
            if self.Name is None:
                name = ''
            else:
                name = str(self.Name)
                if escape_name and not (name.startswith("'") and \
                                        name.endswith("'")):
                    if re.search("""[]['"(),:;_]""", name):
                        name = "'%s'" %  name.replace("'","''")
                    else:
                        name = name.replace(' ','_')
            newick.append(name)
        
        if isinstance(self, PhyloNode):
            if with_distances and self.Length is not None:
                newick.append(":%s" % self.Length)
        
        if semicolon:
            newick.append(";")
        
        return ''.join(newick)

    def getNewick(self, with_distances=False, semicolon=True, escape_name=True):
        """Return the newick string for this tree.

        Arguments:
            - with_distances: whether branch lengths are included.
            - semicolon: end tree string with a semicolon
            - escape_name: if any of these characters []'"(),:;_ exist in a
                nodes name, wrap the name in single quotes

        NOTE: This method returns the Newick representation of this node
        and its descendents. This method is a modification of an implementation
        by Zongzhi Liu
        """
        result = ['(']
        nodes_stack = [[self, len(self.Children)]]
        node_count = 1

        while nodes_stack:
            node_count += 1
            #check the top node, any children left unvisited?
            top = nodes_stack[-1]
            top_node, num_unvisited_children = top
            if num_unvisited_children: #has any child unvisited
                top[1] -= 1  #decrease the #of children unvisited
                next_child = top_node.Children[-num_unvisited_children] # - for order
                #pre-visit
                if next_child.Children:
                    result.append('(')
                nodes_stack.append([next_child, len(next_child.Children)])
            else:  #no unvisited children
                nodes_stack.pop()
                #post-visit
                if top_node.Children:
                    result[-1] = ')'

                if top_node.NameLoaded:
                    if top_node.Name is None:
                        name = ''
                    else:
                        name = str(top_node.Name)
                        if escape_name and not (name.startswith("'") and \
                                                name.endswith("'")):
                            if re.search("""[]['"(),:;_]""", name):
                                name = "'%s'" % name.replace("'", "''")
                            else:
                                name = name.replace(' ','_')
                    result.append(name)
                
                if isinstance(self, PhyloNode):
                    if with_distances and top_node.Length is not None:
                        #result.append(":%s" % top_node.Length)
                        result[-1] = "%s:%s" % (result[-1], top_node.Length)

                result.append(',')

        len_result = len(result)
        if len_result == 2:  # single node no name
            if semicolon:
                return ";"
            else:
                return ''
        elif len_result == 3: # single node with name
            if semicolon:
                return "%s;" % result[1]
            else:
                return result[1]
        else:
            if semicolon:
                result[-1] = ';'
            else:
                result.pop(-1)
            return ''.join(result)
    
    def removeNode(self, target):
        """Removes node by identity instead of value.
        
        Returns True if node was present, False otherwise.
        """
        to_delete = None
        for i, curr_node in enumerate(self.Children):
            if curr_node is target:
                to_delete = i
                break
        if to_delete is None:
            return False
        else:
            del self[to_delete]
            return True
    
    def getEdgeNames(self, tip1name, tip2name,
            getclade, getstem, outgroup_name=None):
        """Return the list of stem and/or sub tree (clade) edge name(s).
        This is done by finding the common intersection, and then getting
        the list of names. If the clade traverses the root, then use the
        outgroup_name argument to ensure valid specification.
        
        Arguments:
            - tip1/2name: edge 1/2 names
            - getstem: whether the name of the clade stem edge is returned.
            - getclade: whether the names of the edges within the clade are
              returned
            - outgroup_name: if provided the calculation is done on a version of
              the tree re-rooted relative to the provided tip.
        
        Usage:
            The returned list can be used to specify subtrees for special
            parameterisation. For instance, say you want to allow the primates
            to have a different value of a particular parameter. In this case,
            provide the results of this method to the parameter controller
            method `setParamRule()` along with the parameter name etc..
        """
        # If outgroup specified put it at the top of the tree so that clades are
        # defined by their distance from it.  This makes a temporary tree with
        # a named edge at it's root, but it's only used here then discarded.
        if outgroup_name is not None:
            outgroup = self.getNodeMatchingName(outgroup_name)
            if outgroup.Children:
                raise TreeError('Outgroup (%s) must be a tip' % outgroup_name)
            self = outgroup.unrootedDeepcopy()
        
        join_edge = self.getConnectingNode(tip1name, tip2name)
        
        edge_names = []
        
        if getstem:
            if join_edge.isroot():
                raise TreeError('LCA(%s,%s) is the root and so has no stem' %
                        (tip1name, tip2name))
            else:
                edge_names.append(join_edge.Name)
        
        if getclade:
            #get the list of names contained by join_edge
            for child in join_edge.Children:
                branchnames = child.getNodeNames(includeself = 1)
                edge_names.extend(branchnames)
        
        return edge_names
    
    def _getNeighboursExcept(self, parent=None):
        # For walking the tree as if it was unrooted.
        return [c for c in (tuple(self.Children) + (self.Parent,))
                if c is not None and c is not parent]
    
    def _getDistances(self, endpoints=None):
        """Iteratively calcluates all of the root-to-tip and tip-to-tip
        distances, resulting in a tuple of:
            - A list of (name, path length) pairs.
            - A dictionary of (tip1,tip2):distance pairs
        """
        ## linearize the tips in postorder.
        # .__start, .__stop compose the slice in tip_order.
        if endpoints is None:
            tip_order = list(self.tips())
        else:
            tip_order = []
            for i,name in enumerate(endpoints):
                node = self.getNodeMatchingName(name)
                tip_order.append(node)
        for i, node in enumerate(tip_order):
            node.__start, node.__stop = i, i+1

        num_tips = len(tip_order)
        result = {}
        tipdistances = zeros((num_tips), float) #distances from tip to curr node

        def update_result():
        # set tip_tip distance between tips of different child
            for child1, child2 in comb(node.Children, 2):
                for tip1 in range(child1.__start, child1.__stop):
                    for tip2 in range(child2.__start, child2.__stop):
                        name1 = tip_order[tip1].Name
                        name2 = tip_order[tip2].Name
                        result[(name1,name2)] = \
                            tipdistances[tip1] + tipdistances[tip2]
                        result[(name2,name1)] = \
                            tipdistances[tip1] + tipdistances[tip2]

        for node in self.traverse(self_before=False, self_after=True):
            if not node.Children:
                continue
            ## subtree with solved child wedges
            starts, stops = [], [] #to calc ._start and ._stop for curr node
            for child in node.Children:
                if hasattr(child, 'Length') and child.Length is not None:
                    child_len = child.Length
                else:
                    child_len = 1 # default length
                tipdistances[child.__start : child.__stop] += child_len
                starts.append(child.__start); stops.append(child.__stop)
            node.__start, node.__stop = min(starts), max(stops)
            ## update result if nessessary
            if len(node.Children) > 1: #not single child
                update_result()

        from_root = []
        for i,n in enumerate(tip_order):
            from_root.append((n.Name, tipdistances[i]))
        return from_root, result

    def getDistances(self, endpoints=None):
        """The distance matrix as a dictionary.
        
        Usage:
            Grabs the branch lengths (evolutionary distances) as
            a complete matrix (i.e. a,b and b,a).
        """
        
        (root_dists, endpoint_dists) = self._getDistances(endpoints)
        return endpoint_dists

    def setMaxTipTipDistance(self):
        """Propagate tip distance information up the tree

        This method was originally implemented by Julia Goodrich with the intent
        of being able to determine max tip to tip distances between nodes on 
        large trees efficiently. The code has been modified to track the 
        specific tips the distance is between
        """
        for n in self.postorder():
            if n.isTip():
                n.MaxDistTips = [[0.0, n.Name], [0.0, n.Name]]
            else:
                if len(n.Children) == 1:
                    tip_a, tip_b = n.Children[0].MaxDistTips
                    tip_a[0] += n.Children[0].Length or 0.0
                    tip_b[0] += n.Children[0].Length or 0.0
                else:
                    tip_info = [(max(c.MaxDistTips), c) for c in n.Children]
                    dists = [i[0][0] for i in tip_info]
                    best_idx = argsort(dists)[-2:]
                    tip_a, child_a = tip_info[best_idx[0]]
                    tip_b, child_b = tip_info[best_idx[1]]
                    tip_a[0] += child_a.Length or 0.0
                    tip_b[0] += child_b.Length or 0.0
                n.MaxDistTips = [tip_a, tip_b]

    def getMaxTipTipDistance(self):
        """Returns the max tip tip distance between any pair of tips
        
        Returns (dist, tip_names, internal_node)
        """
        if not hasattr(self, 'MaxDistTips'):
            self.setMaxTipTipDistance()

        longest = 0.0
        names = [None,None]
        best_node = None
        for n in self.nontips(include_self=True):
            tip_a, tip_b = n.MaxDistTips
            dist = (tip_a[0] + tip_b[0])

            if dist > longest:
                longest = dist
                best_node = n
                names = [tip_a[1], tip_b[1]]
        return longest, names, best_node

    def maxTipTipDistance(self):
        """returns the max distance between any pair of tips
        
        Also returns the tip names  that it is between as a tuple"""
        distmtx, tip_order = self.tipToTipDistances()
        idx_max = divmod(distmtx.argmax(),distmtx.shape[1])
        max_pair = (tip_order[idx_max[0]].Name, tip_order[idx_max[1]].Name)
        return distmtx[idx_max], max_pair
 
    def _getSubTree(self, included_names, constructor=None, keep_root=False):
        """An equivalent node with possibly fewer children, or None"""
        
        # Renumber autonamed edges
        if constructor is None:
            constructor = self._default_tree_constructor()
        
        if self.Name in included_names:
            return self.deepcopy(constructor=constructor)
        else:
            # don't need to pass keep_root to children, though
            # internal nodes will be elminated this way
            children = [child._getSubTree(included_names, constructor)
                    for child in self.Children]
            children = [child for child in children if child is not None]
            if len(children) == 0:
                result = None
            elif len(children) == 1 and not keep_root:
                # Merge parameter dictionaries by adding lengths and making
                # weighted averages of other parameters.  This should probably
                # be moved out of here into a ParameterSet class (Model?) or
                # tree subclass.
                params = {}
                child = children[0]
                if self.Length is not None and child.Length is not None:
                    shared_params = [n for (n,v) in self.params.items()
                        if v is not None
                        and child.params.get(n) is not None
                        and n is not "length"]
                    length = self.Length + child.Length
                    if length:
                        params = dict([(n,
                                (self.params[n]*self.Length +
                                child.params[n]*child.Length) / length)
                            for n in shared_params])
                        params['length'] = length
                result = child
                result.params = params
            else:
                result = constructor(self, tuple(children))
        return result
    
    def getSubTree(self, name_list, ignore_missing=False, keep_root=False):
        """A new instance of a sub tree that contains all the otus that are
        listed in name_list.

        ignore_missing: if False, getSubTree will raise a ValueError if 
        name_list contains names that aren't nodes in the tree

        keep_root: if False, the root of the subtree will be the last common
        ancestor of all nodes kept in the subtree. Root to tip distance is
        then (possibly) different from the original tree
        If True, the root to tip distance remains constant, but root may only
        have one child node.
        """
        edge_names = set(self.getNodeNames(includeself=1, tipsonly=False))
        if not ignore_missing:
            # this may take a long time
            for name in name_list:
                if name not in edge_names:
                    raise ValueError("edge %s not found in tree" % name)
        
        new_tree = self._getSubTree(name_list, keep_root=keep_root)
        if new_tree is None:
            raise TreeError, "no tree created in make sub tree"
        elif new_tree.istip():
            raise TreeError, "only a tip was returned from selecting sub tree"
        else:
            new_tree.Name = "root"
            # keep unrooted
            if len(self.Children) > 2:
                new_tree = new_tree.unrooted()
            return new_tree
    
    def _edgecount(self, parent, cache):
        """"The number of edges beyond 'parent' in the direction of 'self',
        unrooted"""
        neighbours = self._getNeighboursExcept(parent)
        key = (id(parent), id(self))
        if key not in cache:
            cache[key] = 1 + sum([child._edgecount(self, cache)
                    for child in neighbours])
        return cache[key]
    
    def _imbalance(self, parent, cache):
        """The edge count from here, (except via 'parent'), divided into that
        from the heaviest neighbour, and that from the rest of them.  'cache'
        should be a dictionary that can be shared by calls to self.edgecount,
        it stores the edgecount for each node (from self) without having to
        put it on the tree itself."""
        max_weight = 0
        total_weight = 0
        for child in self._getNeighboursExcept(parent):
            weight = child._edgecount(self, cache)
            total_weight += weight
            if weight > max_weight:
                max_weight = weight
                biggest_branch = child
        return (max_weight, total_weight-max_weight, biggest_branch)
    
    def _sorted(self, sort_order):
        """Score all the edges, sort them, and return minimum score and a
        sorted tree.
        """
        # Only need to duplicate whole tree because of .Parent pointers
        
        constructor = self._default_tree_constructor()
        
        if not self.Children:
            tree = self.deepcopy(constructor)
            score = sort_order.index(self.Name)
        else:
            scored_subtrees = [child._sorted(sort_order)
                    for child in self.Children]
            scored_subtrees.sort()
            children = tuple([child.deepcopy(constructor)
                    for (score, child) in scored_subtrees])
            tree = constructor(self, children)
            
            non_null_scores = [score
                    for (score, child) in scored_subtrees if score is not None]
            score = (non_null_scores or [None])[0]
        return (score, tree)
    
    def sorted(self, sort_order=[]):
        """An equivalent tree sorted into a standard order. If this is not
        specified then alphabetical order is used.  At each node starting from
        root, the algorithm will try to put the descendant which contains the
        lowest scoring tip on the left.
        """
        tip_names = self.getTipNames()
        tip_names.sort()
        full_sort_order = sort_order + tip_names
        (score, tree) = self._sorted(full_sort_order)
        return tree
    
    def _asciiArt(self, char1='-', show_internal=True, compact=False):
        LEN = 10
        PAD = ' ' * LEN
        PA = ' ' * (LEN-1)
        namestr = self.Name or '' # prevents name of NoneType
        if self.Children:
            mids = []
            result = []
            for c in self.Children:
                if c is self.Children[0]:
                    char2 = '/'
                elif c is self.Children[-1]:
                    char2 = '\\'
                else:
                    char2 = '-'
                (clines, mid) = c._asciiArt(char2, show_internal, compact)
                mids.append(mid+len(result))
                result.extend(clines)
                if not compact:
                    result.append('')
            if not compact:
                result.pop()
            (lo, hi, end) = (mids[0], mids[-1], len(result))
            prefixes = [PAD] * (lo+1) + [PA+'|'] * (hi-lo-1) + [PAD] * (end-hi)
            mid = (lo + hi) / 2
            prefixes[mid] = char1 + '-'*(LEN-2) + prefixes[mid][-1]
            result = [p+l for (p,l) in zip(prefixes, result)]
            if show_internal:
                stem = result[mid]
                result[mid] = stem[0] + namestr + stem[len(namestr)+1:]
            return (result, mid)
        else:
            return ([char1 + '-' + namestr], 0)
    
    def asciiArt(self, show_internal=True, compact=False):
        """Returns a string containing an ascii drawing of the tree.
        
        Arguments:
        - show_internal: includes internal edge names.
        - compact: use exactly one line per tip.
        """
        (lines, mid) = self._asciiArt(
                show_internal=show_internal, compact=compact)
        return '\n'.join(lines)
    
    def _getXmlLines(self, indent=0, parent_params=None):
        """Return the xml strings for this edge.
        """
        params = {}
        if parent_params is not None:
            params.update(parent_params)
        pad = '  ' * indent
        xml = ["%s<clade>" % pad]
        if self.NameLoaded:
            xml.append("%s   <name>%s</name>" % (pad, self.Name))
        for (n,v) in self.params.items():
            if v == params.get(n, None):
                continue
            xml.append("%s   <param><name>%s</name><value>%s</value></param>"
                    % (pad, n, v))
            params[n] = v
        for child in self.Children:
            xml.extend(child._getXmlLines(indent + 1, params))
        xml.append(pad + "</clade>")
        return xml
    
    def getXML(self):
        """Return XML formatted tree string."""
        header = ['<?xml version="1.0"?>']  # <!DOCTYPE ...
        return '\n'.join(header + self._getXmlLines())
    
    def writeToFile(self, filename, with_distances=True, format=None):
        """Save the tree to filename
        
        Arguments:
            - filename: self-evident
            - with_distances: whether branch lengths are included in string.
            - format: default is newick, xml is alternate. Argument overrides
              the filename suffix. All attributes are saved in the xml format.
        """
        if format:
            xml = format.lower() == 'xml'
        else:
            xml = filename.lower().endswith('xml')
        
        if xml:
            data = self.getXML()
        else:
            data = self.getNewick(with_distances=with_distances)
        outf = open(filename, "w")
        outf.writelines(data)
        outf.close()

    def getNodeNames(self, includeself=True, tipsonly=False):
        """Return a list of edges from this edge - may or may not include self.
        This node (or first connection) will be the first, and then they will
        be listed in the natural traverse order.
        """
        if tipsonly:
            nodes = self.traverse(self_before=False, self_after=False)
        else:
            nodes = list(self.traverse())
            if not includeself:
                nodes = nodes[:-1]
        return [node.Name for node in nodes]
   
    def getTipNames(self, includeself=False):
        """return the list of the names of all tips contained by this edge
        """
        return self.getNodeNames(includeself, tipsonly=True)
    
    def getEdgeVector(self):
        """Collect the list of edges in postfix order"""
        return [node for node in self.traverse(False, True)]
    
    def _getNodeMatchingName(self, name):
        """
        find the edge with the name, or return None
        """
        for node in self.traverse(self_before=True, self_after=False):
            if node.Name == name:
                return node
        return None
    
    def getNodeMatchingName(self, name):
        node = self._getNodeMatchingName(name)
        if node is None:
            raise TreeError("No node named '%s' in %s" %
                    (name, self.getTipNames()))
        return node
   
    def getConnectingNode(self, name1, name2):
        """Finds the last common ancestor of the two named edges."""
        edge1 = self.getNodeMatchingName(name1)
        edge2 = self.getNodeMatchingName(name2)
        lca = edge1.lastCommonAncestor(edge2)
        if lca is None:
            raise TreeError("No LCA found for %s and %s" % (name1, name2))
        return lca
    
    def getConnectingEdges(self, name1, name2):
        """returns a list of edges connecting two nodes
    
        includes self and other in the list"""
        edge1 = self.getNodeMatchingName(name1)
        edge2 = self.getNodeMatchingName(name2)
        LCA = self.getConnectingNode(name1, name2)
        node_path = [edge1]
        node_path.extend(edge1.ancestors())
        #remove nodes deeper than the LCA
        LCA_ind = node_path.index(LCA)
        node_path = node_path[:LCA_ind+1]
        #remove LCA and deeper nodes from anc list of other
        anc2 = edge2.ancestors()
        LCA_ind = anc2.index(LCA)
        anc2 = anc2[:LCA_ind]
        anc2.reverse()
        node_path.extend(anc2)
        node_path.append(edge2)
        return node_path
    
    def getParamValue(self, param, edge):
        """returns the parameter value for named edge"""
        return self.getNodeMatchingName(edge).params[param]
    
    def setParamValue(self, param, edge, value):
        """set's the value for param at named edge"""
        self.getNodeMatchingName(edge).params[param] = value

    def reassignNames(self, mapping, nodes=None):
        """Reassigns node names based on a mapping dict

        mapping : dict, old_name -> new_name
        nodes : specific nodes for renaming (such as just tips, etc...)
        """
        if nodes is None:
            nodes = self.traverse()

        for n in nodes:
            if n.Name in mapping:
                n.Name = mapping[n.Name]

    def getNodesDict(self):
        """Returns a dict keyed by node name, value is node

        Will raise TreeError if non-unique names are encountered
        """
        res = {}

        for n in self.traverse():
            if n.Name in res:
                raise TreeError, "getNodesDict requires unique node names"
            else:
                res[n.Name] = n

        return res

    def subset(self):
        """Returns set of names that descend from specified node"""
        return frozenset([i.Name for i in self.tips()])

    def subsets(self):
        """Returns all sets of names that come from specified node and its kids"""
        sets = []
        for i in self.traverse(self_before=False, self_after=True, \
            include_self=False):
            if not i.Children:
                i.__leaf_set = frozenset([i.Name])
            else:
                leaf_set = reduce(or_, [c.__leaf_set for c in i.Children])
                if len(leaf_set) > 1:
                    sets.append(leaf_set)
                i.__leaf_set = leaf_set
        return frozenset(sets)
                
    def compareBySubsets(self, other, exclude_absent_taxa=False):
        """Returns fraction of overlapping subsets where self and other differ.
            
        Other is expected to be a tree object compatible with PhyloNode.
        
        Note: names present in only one of the two trees will count as 
        mismatches: if you don't want this behavior, strip out the non-matching
        tips first.
        """
        self_sets, other_sets = self.subsets(), other.subsets()
        if exclude_absent_taxa:
            in_both = self.subset() & other.subset()
            self_sets = [i & in_both for i in self_sets]
            self_sets = frozenset([i for i in self_sets if len(i) > 1])
            other_sets = [i & in_both for i in other_sets]
            other_sets = frozenset([i for i in other_sets if len(i) > 1])
        total_subsets = len(self_sets) + len(other_sets)
        intersection_length = len(self_sets & other_sets)
        if not total_subsets:   #no common subsets after filtering, so max dist
            return 1
        return 1 - 2*intersection_length/float(total_subsets)

    def tipToTipDistances(self, default_length=1):
        """Returns distance matrix between all pairs of tips, and a tip order.
            
        Warning: .__start and .__stop added to self and its descendants.

        tip_order contains the actual node objects, not their names (may be
        confusing in some cases).
        """
        ## linearize the tips in postorder.
        # .__start, .__stop compose the slice in tip_order.
        tip_order = list(self.tips())
        for i, tip in enumerate(tip_order):
            tip.__start, tip.__stop = i, i+1

        num_tips = len(tip_order)
        result = zeros((num_tips, num_tips), float) #tip by tip matrix
        tipdistances = zeros((num_tips), float) #distances from tip to curr node

        def update_result(): 
        # set tip_tip distance between tips of different child
            for child1, child2 in comb(node.Children, 2):
                for tip1 in range(child1.__start, child1.__stop):
                    for tip2 in range(child2.__start, child2.__stop):
                        result[tip1,tip2] = \
                            tipdistances[tip1] + tipdistances[tip2]

        for node in self.traverse(self_before=False, self_after=True):
            if not node.Children:
                continue
            ## subtree with solved child wedges
            starts, stops = [], [] #to calc ._start and ._stop for curr node
            for child in node.Children:
                if hasattr(child, 'Length') and child.Length is not None:
                    child_len = child.Length
                else:
                    child_len = default_length
                tipdistances[child.__start : child.__stop] += child_len
                starts.append(child.__start); stops.append(child.__stop)
            node.__start, node.__stop = min(starts), max(stops)
            ## update result if nessessary
            if len(node.Children) > 1: #not single child
                update_result()
        return result+result.T, tip_order 

    def compareByTipDistances(self, other, dist_f=distance_from_r):
        """Compares self to other using tip-to-tip distance matrices.

        Value returned is dist_f(m1, m2) for the two matrices. Default is
        to use the Pearson correlation coefficient, with +1 giving a distance
        of 0 and -1 giving a distance of +1 (the madimum possible value).
        Depending on the application, you might instead want to use
        distance_from_r_squared, which counts correlations of both +1 and -1
        as identical (0 distance).
        
        Note: automatically strips out the names that don't match (this is
        necessary for this method because the distance between non-matching 
        names and matching names is undefined in the tree where they don't 
        match, and because we need to reorder the names in the two trees to 
        match up the distance matrices).
        """
        self_names = [i.Name for i in self.tips()]
        other_names = [i.Name for i in other.tips()]
        common_names = frozenset(self_names) & frozenset(other_names)
        if not common_names:
            raise ValueError, "No names in common between the two trees."""
        if len(common_names) <= 2:
            return 1    #the two trees must match by definition in this case
        #figure out correct order of the two name matrices
        self_order = [self_names.index(i) for i in common_names]
        other_order = [other_names.index(i) for i in common_names]
        self_matrix = self.tipToTipDistances()[0][self_order][:,self_order]
        other_matrix = other.tipToTipDistances()[0][other_order][:,other_order]
        return dist_f(self_matrix, other_matrix)

class PhyloNode(TreeNode):

    def __init__(self, *args, **kwargs):
        length = kwargs.get('Length', None)
        params = kwargs.get('Params', {})
        if 'length' not in params:
            params['length'] = length
        kwargs['Params'] = params
        super(PhyloNode, self).__init__(*args, **kwargs)

    def _set_length(self, value):
        if not hasattr(self, "params"):
            self.params = {}
        self.params["length"] = value
    
    def _get_length(self):
        return self.params.get("length", None)
    
    Length = property(_get_length, _set_length)

    def getNewick(self, with_distances=False, semicolon=True, escape_name=True):
        return TreeNode.getNewick(self, with_distances, semicolon, escape_name)
    
    def __str__(self):
        """Returns string version of self, with names and distances."""
        return self.getNewick(with_distances=True)
    
    
    def distance(self, other):
        """Returns branch length between self and other."""
        #never any length between self and other
        if self is other:
            return 0
        #otherwise, find self's ancestors and find the first ancestor of
        #other that is in the list
        self_anc = self.ancestors()
        self_anc_dict = dict([(id(n),n) for n in self_anc])
        self_anc_dict[id(self)] = self
        
        count = 0
        while other is not None:
            if id(other) in self_anc_dict:
                #found the first shared ancestor -- need to sum other branch
                curr = self
                while curr is not other:
                    if curr.Length:
                        count += curr.Length
                    curr = curr._parent
                return count
            else:
                if other.Length:
                    count += other.Length
            other = other._parent
        return None

    def totalDescendingBranchLength(self):
        """Returns total descending branch length from self"""
        return sum([n.Length for n in self.traverse(include_self=False) \
                                     if n.Length is not None])

    def tipsWithinDistance(self, distance):
        """Returns tips within specified distance from self

        Branch lengths of None will be interpreted as 0
        """
        def get_distance(d1, d2):
            if d2 is None:
                return d1
            else:
                return d1 + d2

        to_process = [(self, 0.0)]
        tips_to_save = []

        curr_node, curr_dist = to_process[0]

        seen = set([id(self)])
        while to_process:
            curr_node, curr_dist = to_process.pop(0)
           
            # have we've found a tip within distance?
            if curr_node.isTip() and curr_node != self:
                tips_to_save.append(curr_node)
                continue
        
            # add the parent node if it is within distance
            parent_dist = get_distance(curr_dist, curr_node.Length)
            if curr_node.Parent is not None and parent_dist <= distance and \
                    id(curr_node.Parent) not in seen:
                to_process.append((curr_node.Parent, parent_dist))
                seen.add(id(curr_node.Parent))

            # add children if we haven't seen them and if they are in distance
            for child in curr_node.Children:
                if id(child) in seen:
                    continue
                seen.add(id(child))

                child_dist = get_distance(curr_dist, child.Length)
                if child_dist <= distance:
                    to_process.append((child, child_dist))

        return tips_to_save

    def prune(self):
        """Reconstructs correct tree after nodes have been removed.
        
        Internal nodes with only one child will be removed and new connections
        and Branch lengths will be made to reflect change. 
        """
        #traverse tree to decide nodes to be removed.
        nodes_to_remove = []
        for node in self.traverse():
            if (node.Parent is not None) and (len(node.Children)==1):
                nodes_to_remove.append(node)
        for node in nodes_to_remove:
            #save current parent
            curr_parent=node.Parent
            #save child
            child=node.Children[0]
            #remove current node by setting parent to None
            node.Parent=None
            #Connect child to current node's parent
            child.Parent=curr_parent
            #Add the Length of the removed node to the Length of the Child
            if child.Length is None or node.Length is None:
                child.Length = child.Length or node.Length
            else:
                child.Length = child.Length + node.Length

    def unrootedDeepcopy(self, constructor=None, parent=None):
        # walks the tree unrooted-style, ie: treating self.Parent as just
        # another child 'parent' is where we got here from, ie: the neighbour
        # that we don't need to explore.
        if constructor is None:
            constructor = self._default_tree_constructor()
        
        neighbours = self._getNeighboursExcept(parent)
        children = []
        for child in neighbours:
            children.append(child.unrootedDeepcopy(constructor, parent=self))
        
        # we might be walking UP the tree, so:
        if parent is None:
            # base edge
            edge = None
        elif parent.Parent is self:
            # self's parent is becoming self's child, and edge params are stored
            # by the child
            edge = parent
        else:
            assert parent is self.Parent
            edge = self
        
        result = constructor(edge, tuple(children))
        if parent is None:
            result.Name = "root"
        return result
    
    def bifurcating(self, constructor=None):
        # With every node having 2 or fewer children.
        if constructor is None:
            constructor = self._default_tree_constructor()
        children = [child.bifurcating(constructor) for child in self.Children]
        while len(children) > 2:
            extra = constructor(None, tuple(children[-2:]))
            children[-2:] = [extra]
        result = constructor(self, tuple(children))
        return result
    
    def balanced(self):
        """Tree 'rooted' here with no neighbour having > 50% of the edges.
        
        Usage:
            Using a balanced tree can substantially improve performance of
            the likelihood calculations. Note that the resulting tree has a
            different orientation with the effect that specifying clades or
            stems for model parameterisation should be done using the
            'outgroup_name' argument.
        """
        # this should work OK on ordinary 3-way trees, not so sure about
        # other cases.  Given 3 neighbours, if one has > 50% of edges it
        # can only improve things to divide it up, worst case:
        # (51),25,24 -> (50,1),49.
        # If no neighbour has >50% we can't improve on where we are, eg:
        # (49),25,26 -> (20,19),51
        last_edge = None
        edge = self
        known_weight = 0
        cache = {}
        while 1:
            (max_weight, remaining_weight, next_edge) = edge._imbalance(
                    last_edge, cache)
            known_weight += remaining_weight
            if max_weight <= known_weight+2:
                break
            last_edge = edge
            edge = next_edge
            known_weight += 1
        return edge.unrootedDeepcopy()
    
    def sameTopology(self, other):
        """Tests whether two trees have the same topology."""
        tip_names = self.getTipNames()
        root_at = tip_names[0]
        me = self.rootedWithTip(root_at).sorted(tip_names)
        them = other.rootedWithTip(root_at).sorted(tip_names)
        return self is other or me.sameShape(them)
    
    def unrooted(self):
        """A tree with at least 3 children at the root.
        """
        constructor = self._default_tree_constructor()
        need_to_expand = len(self.Children) < 3
        new_children = []
        for oldnode in self.Children:
            if oldnode.Children and need_to_expand:
                for sib in oldnode.Children:
                    sib = sib.deepcopy(constructor)
                    if sib.Length is not None and oldnode.Length is not None:
                        sib.Length += oldnode.Length
                    new_children.append(sib)
                need_to_expand = False
            else:
                new_children.append(oldnode.deepcopy(constructor))
        return constructor(self, new_children)
    
    def rootedAt(self, edge_name):
        """Return a new tree rooted at the provided node.
        
        Usage:
            This can be useful for drawing unrooted trees with an orientation
            that reflects knowledge of the true root location.
        """
        newroot = self.getNodeMatchingName(edge_name)
        if not newroot.Children:
            raise TreeError("Can't use a tip (%s) as the root" %
                    repr(edge_name))
        return newroot.unrootedDeepcopy()
    
    def rootedWithTip(self, outgroup_name):
        """A new tree with the named tip as one of the root's children"""
        tip = self.getNodeMatchingName(outgroup_name)
        return tip.Parent.unrootedDeepcopy()

    def rootAtMidpoint(self):
        """ return a new tree rooted at midpoint of the two tips farthest apart
    
        this fn doesn't preserve the internal node naming or structure,
        but does keep tip to tip distances correct.  uses unrootedDeepcopy()
        """
        # max_dist, tip_names = tree.maxTipTipDistance()
        # this is slow


        max_dist, tip_names = self.maxTipTipDistance()
        half_max_dist = max_dist/2.0
        if max_dist == 0.0: # only pathological cases with no lengths
            return self.unrootedDeepcopy()
        # print tip_names
        tip1 = self.getNodeMatchingName(tip_names[0])
        tip2 = self.getNodeMatchingName(tip_names[1])
        lca = self.getConnectingNode(tip_names[0],tip_names[1]) # last comm ancestor
        if tip1.distance(lca) > half_max_dist:
            climb_node = tip1
        else:
            climb_node = tip2
        
        dist_climbed = 0.0
        while dist_climbed + climb_node.Length < half_max_dist:
            dist_climbed += climb_node.Length
            climb_node = climb_node.Parent
            
        # now midpt is either at on the branch to climb_node's  parent
        # or midpt is at climb_node's parent
        # print dist_climbed, half_max_dist, 'dists cl hamax'
        if dist_climbed + climb_node.Length == half_max_dist:
            # climb to midpoint spot
            climb_node = climb_node.Parent
            if climb_node.isTip():
                raise RuntimeError('error trying to root tree at tip')
            else:
                # print climb_node.Name, 'clmb node'
                return climb_node.unrootedDeepcopy()
        
        else:
            # make a new node on climb_node's branch to its parent
            old_br_len = climb_node.Length
            new_root = type(self)()
            new_root.Parent = climb_node.Parent
            climb_node.Parent = new_root
            climb_node.Length = half_max_dist - dist_climbed
            new_root.Length = old_br_len - climb_node.Length
            return new_root.unrootedDeepcopy()

    
    def _find_midpoint_nodes(self, max_dist, tip_pair):
        """returns the nodes surrounding the maxTipTipDistance midpoint 
        
        WAS used for midpoint rooting.  ORPHANED NOW
        max_dist: The maximum distance between any 2 tips
        tip_pair: Names of the two tips associated with max_dist
        """
        half_max_dist = max_dist/2.0
        #get a list of the nodes that separate the tip pair
        node_path = self.getConnectingEdges(tip_pair[0], tip_pair[1])
        tip1 = self.getNodeMatchingName(tip_pair[0])
        for index, node in enumerate(node_path):
            dist = tip1.distance(node)
            if dist > half_max_dist:
                return node, node_path[index-1]

    def setTipDistances(self):
        """Sets distance from each node to the most distant tip."""
        for node in self.traverse(self_before=False, self_after=True):
            if node.Children:
                node.TipDistance = max([c.Length + c.TipDistance for \
                    c in node.Children])
            else:
                node.TipDistance = 0

    def scaleBranchLengths(self, max_length=100, ultrametric=False):
        """Scales BranchLengths in place to integers for ascii output.

        Warning: tree might not be exactly the length you specify.

        Set ultrametric=True if you want all the root-tip distances to end
        up precisely the same.
        """
        self.setTipDistances()
        orig_max = max([n.TipDistance for n in self.traverse()])
        if not ultrametric: #easy case -- just scale and round
            for node in self.traverse():
                curr = node.Length
                if curr is not None:
                    node.ScaledBranchLength =  \
                        max(1, int(round(1.0*curr/orig_max*max_length)))
        else:   #hard case -- need to make sure they all line up at the end
            for node in self.traverse(self_before=False, self_after=True):
                if not node.Children:   #easy case: ignore tips
                    node.DistanceUsed = 0
                    continue
                #if we get here, we know the node has children
                #figure out what distance we want to set for this node
                ideal_distance=int(round(node.TipDistance/orig_max*max_length))
                min_distance = max([c.DistanceUsed for c in node.Children]) + 1
                distance = max(min_distance, ideal_distance)
                for c in node.Children:
                    c.ScaledBranchLength = distance - c.DistanceUsed
                node.DistanceUsed = distance
        #reset the BranchLengths
        for node in self.traverse(self_before=True, self_after=False):
            if node.Length is not None:
                node.Length = node.ScaledBranchLength
            if hasattr(node, 'ScaledBranchLength'):
                del node.ScaledBranchLength
            if hasattr(node, 'DistanceUsed'):
                del node.DistanceUsed
            if hasattr(node, 'TipDistance'):
                del node.TipDistance

    def _getDistances(self, endpoints=None):
        """Iteratively calcluates all of the root-to-tip and tip-to-tip
        distances, resulting in a tuple of:
            - A list of (name, path length) pairs.
            - A dictionary of (tip1,tip2):distance pairs
        """
        ## linearize the tips in postorder.
        # .__start, .__stop compose the slice in tip_order.
        if endpoints is None:
            tip_order = list(self.tips())
        else:
            tip_order = []
            for i,name in enumerate(endpoints):
                node = self.getNodeMatchingName(name)
                tip_order.append(node)
        for i, node in enumerate(tip_order):
            node.__start, node.__stop = i, i+1

        num_tips = len(tip_order)
        result = {}
        tipdistances = zeros((num_tips), float) #distances from tip to curr node

        def update_result():
        # set tip_tip distance between tips of different child
            for child1, child2 in comb(node.Children, 2):
                for tip1 in range(child1.__start, child1.__stop):
                    for tip2 in range(child2.__start, child2.__stop):
                        name1 = tip_order[tip1].Name
                        name2 = tip_order[tip2].Name
                        result[(name1,name2)] = \
                            tipdistances[tip1] + tipdistances[tip2]
                        result[(name2,name1)] = \
                            tipdistances[tip1] + tipdistances[tip2]

        for node in self.traverse(self_before=False, self_after=True):
            if not node.Children:
                continue
            ## subtree with solved child wedges
            starts, stops = [], [] #to calc ._start and ._stop for curr node
            for child in node.Children:
                if hasattr(child, 'Length') and child.Length is not None:
                    child_len = child.Length
                else:
                    child_len = 1 # default length
                tipdistances[child.__start : child.__stop] += child_len
                starts.append(child.__start); stops.append(child.__stop)
            node.__start, node.__stop = min(starts), max(stops)
            ## update result if nessessary
            if len(node.Children) > 1: #not single child
                update_result()

        from_root = []
        for i,n in enumerate(tip_order):
            from_root.append((n.Name, tipdistances[i]))
        return from_root, result

    def getDistances(self, endpoints=None):
        """The distance matrix as a dictionary.
        
        Usage:
            Grabs the branch lengths (evolutionary distances) as
            a complete matrix (i.e. a,b and b,a)."""

        (root_dists, endpoint_dists) = self._getDistances(endpoints)
        return endpoint_dists

    def tipToTipDistances(self, endpoints=None, default_length=1):
        """Returns distance matrix between all pairs of tips, and a tip order.
            
        Warning: .__start and .__stop added to self and its descendants.

        tip_order contains the actual node objects, not their names (may be
        confusing in some cases).
        """
        all_tips = self.tips()
        if endpoints is None:
            tip_order = list(all_tips)
        else:
             if isinstance(endpoints[0], PhyloNode):
                tip_order = endpoints
             else:
                tip_order = [self.getNodeMatchingName(n) for n in endpoints]

        ## linearize all tips in postorder
        # .__start, .__stop compose the slice in tip_order.
        for i, node in enumerate(all_tips):
            node.__start, node.__stop = i, i+1
        
        # the result map provides index in the result matrix
        result_map = dict([(n.__start,i) for i,n in enumerate(tip_order)])
        num_all_tips = len(all_tips) # total number of tips
        num_tips = len(tip_order) # total number of tips in result
        result = zeros((num_tips, num_tips), float) # tip by tip matrix
        tipdistances = zeros((num_all_tips), float) # dist from tip to curr node

        def update_result():
        # set tip_tip distance between tips of different child
            for child1, child2 in comb(node.Children, 2):
                for tip1 in range(child1.__start, child1.__stop):
                    if tip1 not in result_map:
                        continue
                    res_tip1 = result_map[tip1]
                    for tip2 in range(child2.__start, child2.__stop):
                        if tip2 not in result_map:
                            continue
                        result[res_tip1,result_map[tip2]] = \
                            tipdistances[tip1] + tipdistances[tip2]

        for node in self.traverse(self_before=False, self_after=True):
            if not node.Children:
                continue
            ## subtree with solved child wedges
            starts, stops = [], [] #to calc ._start and ._stop for curr node
            for child in node.Children:
                if hasattr(child, 'Length') and child.Length is not None:
                    child_len = child.Length
                else:
                    child_len = default_length
                tipdistances[child.__start : child.__stop] += child_len
                starts.append(child.__start); stops.append(child.__stop)
            node.__start, node.__stop = min(starts), max(stops)
            ## update result if nessessary
            if len(node.Children) > 1: #not single child
                update_result()
        return result+result.T, tip_order

    def compareByTipDistances(self, other, sample=None, dist_f=distance_from_r,\
            shuffle_f=shuffle):
        """Compares self to other using tip-to-tip distance matrices.

        Value returned is dist_f(m1, m2) for the two matrices. Default is
        to use the Pearson correlation coefficient, with +1 giving a distance
        of 0 and -1 giving a distance of +1 (the madimum possible value).
        Depending on the application, you might instead want to use
        distance_from_r_squared, which counts correlations of both +1 and -1
        as identical (0 distance).
        
        Note: automatically strips out the names that don't match (this is
        necessary for this method because the distance between non-matching 
        names and matching names is undefined in the tree where they don't 
        match, and because we need to reorder the names in the two trees to 
        match up the distance matrices).
        """
        self_names = dict([(i.Name, i) for i in self.tips()])
        other_names = dict([(i.Name, i) for i in other.tips()])
        common_names = frozenset(self_names.keys()) & \
                       frozenset(other_names.keys())
        common_names = list(common_names)

        if not common_names:
            raise ValueError, "No names in common between the two trees."""
        if len(common_names) <= 2:
            return 1    #the two trees must match by definition in this case

        if sample is not None:
            shuffle_f(common_names)
            common_names = common_names[:sample]
            
        self_nodes = [self_names[k] for k in common_names]
        other_nodes = [other_names[k] for k in common_names]

        self_matrix = self.tipToTipDistances(endpoints=self_nodes)[0]
        other_matrix = other.tipToTipDistances(endpoints=other_nodes)[0]

        return dist_f(self_matrix, other_matrix)

class TreeBuilder(object):
    # Some tree code which isn't needed once the tree is finished.
    # Mostly exists to give edges unique names
    # Children must be created before their parents.
    
    def __init__(self, mutable=False, constructor=PhyloNode):
        self._used_names = {'edge':-1}
        self._known_edges = {}
        self.TreeNodeClass = constructor
    
    def _unique_name(self, name):
        # Unnamed edges become edge.0, edge.1 edge.2 ...
        # Other duplicates go mouse mouse.2 mouse.3 ...
        if not name:
            name = 'edge'
        if name in self._used_names:
            self._used_names[name] += 1
            name += '.' + str(self._used_names[name])
            name = self._unique_name(name) # in case of names like 'edge.1.1'
        else:
            self._used_names[name] = 1
        return name
    
    def _params_for_edge(self, edge):
        # default is just to keep it
        return edge.params
    
    def edgeFromEdge(self, edge, children, params=None):
        """Callback for tree-to-tree transforms like getSubTree"""
        if edge is None:
            assert not params
            return self.createEdge(children, "root", {}, False)
        else:
            if params is None:
                params = self._params_for_edge(edge)
            return self.createEdge(
                    children, edge.Name, params, nameLoaded=edge.NameLoaded)
    
    def createEdge(self, children, name, params, nameLoaded=True):
        """Callback for newick parser"""
        if children is None:
            children = []
        node = self.TreeNodeClass(
                Children = list(children),
                Name = self._unique_name(name),
                NameLoaded = nameLoaded and (name is not None),
                Params = params,
                )
        self._known_edges[id(node)] = node
        return node
    

### from qiime, but using this tree as default ctor
from cogent.parse.tree import DndParser
def parse_newick(lines, constructor=PhyloNode):
    """Return PhyloNode from newick file handle stripping quotes from tip names
    
        This function wraps cogent.parse.tree.DndParser stripping 
         matched leading/trailing single quotes from tip names, and returning
         a PhyloNode object by default (alternate constructor can be passed 
         with constructor=).
         
        Sripping of quotes is essential for many applications in Qiime, as 
         the tip names are frequently matched to OTU ids, and if the tip name
         is read in with leading/trailing quotes, node.Name won't match to the
         corresponding OTU identifier. Disaster follows.
        
    """
    return DndParser(lines, constructor=constructor, unescape_name=True)