# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

# This file was automatically generated from ply. To re-generate this file,
# remove it from this folder, then build astropy and run the tests in-place:
#
#   python setup.py build_ext --inplace
#   pytest astroquery/simbad
#
# You can then commit the changes to this file.


# criteria_parsetab.py
# This file is automatically generated. Do not edit.
# pylint: disable=W,C,R
_tabversion = '3.10'

_lr_method = 'LALR'

_lr_signature = "BINARY_OPERATOR COLUMN IN LIKE LIST NOTLIKE NUMBER REGION STRINGcriteria : criteria '|' criteriacriteria : criteria '&' criteriacriteria : '(' criteria ')'criteria : COLUMN BINARY_OPERATOR STRING\n                        | COLUMN BINARY_OPERATOR NUMBER\n                        | COLUMN IN LIST\n            criteria : COLUMN BINARY_OPERATOR COLUMN\n            criteria : COLUMN LIKE STRINGcriteria : COLUMN NOTLIKE STRINGcriteria : REGION"
    
_lr_action_items = {'(':([0,2,5,6,],[2,2,2,2,]),'COLUMN':([0,2,5,6,8,],[3,3,3,3,15,]),'REGION':([0,2,5,6,],[4,4,4,4,]),'$end':([1,4,12,13,14,15,16,17,18,19,20,],[0,-10,-1,-2,-3,-7,-4,-5,-6,-8,-9,]),'|':([1,4,7,12,13,14,15,16,17,18,19,20,],[5,-10,5,5,5,-3,-7,-4,-5,-6,-8,-9,]),'&':([1,4,7,12,13,14,15,16,17,18,19,20,],[6,-10,6,6,6,-3,-7,-4,-5,-6,-8,-9,]),'BINARY_OPERATOR':([3,],[8,]),'IN':([3,],[9,]),'LIKE':([3,],[10,]),'NOTLIKE':([3,],[11,]),')':([4,7,12,13,14,15,16,17,18,19,20,],[-10,14,-1,-2,-3,-7,-4,-5,-6,-8,-9,]),'STRING':([8,10,11,],[16,19,20,]),'NUMBER':([8,],[17,]),'LIST':([9,],[18,]),}

_lr_action = {}
for _k, _v in _lr_action_items.items():
   for _x,_y in zip(_v[0],_v[1]):
      if not _x in _lr_action:  _lr_action[_x] = {}
      _lr_action[_x][_k] = _y
del _lr_action_items

_lr_goto_items = {'criteria':([0,2,5,6,],[1,7,12,13,]),}

_lr_goto = {}
for _k, _v in _lr_goto_items.items():
   for _x, _y in zip(_v[0], _v[1]):
       if not _x in _lr_goto: _lr_goto[_x] = {}
       _lr_goto[_x][_k] = _y
del _lr_goto_items
_lr_productions = [
  ("S' -> criteria","S'",1,None,None,None),
  ('criteria -> criteria | criteria','criteria',3,'p_criteria_OR','utils.py',298),
  ('criteria -> criteria & criteria','criteria',3,'p_criteria_AND','utils.py',302),
  ('criteria -> ( criteria )','criteria',3,'p_criteria_parenthesis','utils.py',306),
  ('criteria -> COLUMN BINARY_OPERATOR STRING','criteria',3,'p_criteria_string','utils.py',310),
  ('criteria -> COLUMN BINARY_OPERATOR NUMBER','criteria',3,'p_criteria_string','utils.py',311),
  ('criteria -> COLUMN IN LIST','criteria',3,'p_criteria_string','utils.py',312),
  ('criteria -> COLUMN BINARY_OPERATOR COLUMN','criteria',3,'p_criteria_string_no_ticks','utils.py',317),
  ('criteria -> COLUMN LIKE STRING','criteria',3,'p_criteria_like','utils.py',323),
  ('criteria -> COLUMN NOTLIKE STRING','criteria',3,'p_criteria_notlike','utils.py',327),
  ('criteria -> REGION','criteria',1,'p_criteria_region','utils.py',331),
]
