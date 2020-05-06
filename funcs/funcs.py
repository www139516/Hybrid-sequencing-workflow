#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time : 2020/5/6 9:41
# @Author : Yuancong Wang
# @Site : 
# @File : funcs.py
# @Software: PyCharm
'''
Description: This python contains common functions used by different classes
Output: Depends
Input: Depends
Other notes: None
'''


import os


def get_abspath_from_input_path(in_path):

    wkdir = os.getcwd()
    out_path = ''
    if not in_path:
        out_path = wkdir
    elif '/' not in in_path:
        out_path = os.path.join(wkdir, in_path)
    else:
        out_path = in_path
    return os.path.abspath(out_path)

