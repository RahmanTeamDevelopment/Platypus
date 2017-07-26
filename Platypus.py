#!/usr/bin/python

from __future__ import division

import sys
import platypus.runner


if __name__ == "__main__":

    commands = {}
    commands["callVariants"] = platypus.runner.callVariants
    commands["printRegions"] = platypus.runner.printRegions

    if len(sys.argv) == 1 or sys.argv[1] not in commands.keys():
        print "\n\n"
        print "Invalid usage: use Platypus as follows:"
        print "\n"
        for k in commands.keys():
            print "python Platypus.py", k, "[Options]"
        print "\n"
        print "For a list of possible options for a specific command, type 'python Platypus.py Command -h'"
        print "\n"
        sys.exit(0)
    if sys.version_info[0] != 2 or sys.version_info[1] < 6:
        print "\n\n"
        print "Platypus works only with Python versions 2.6 and greater. Python 3.X is not yet supported."
        print "\n\n"
        sys.exit(0)
    else:
        command = sys.argv[1]
        commands[command](sys.argv[2:])
