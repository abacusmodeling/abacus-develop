'''
Date: 2021-08-21 11:43:27
LastEditors: jiyuyang
LastEditTime: 2021-08-21 22:03:06
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

import argparse

from abacus_plot.cmdline import Show


def main():
    parser = argparse.ArgumentParser(
        prog='abacus-plot', description='Plotting tools for ABACUS')

    # Show
    parser.add_argument('-b', '--band', dest='band', type=str,
                        default=None, help='plot band structure and show band information.')
    parser.add_argument('-d', '--dos', dest='dos', type=str,
                        default=None, help='plot density of state(DOS).')
    parser.set_defaults(func=Show().show_cmdline)

    args = parser.parse_args()
    args.func(args)
