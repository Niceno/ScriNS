"""
Run all test
"""

import types

from scrins import tests

def main():
    test_mods = [obj for name, obj in vars(tests).items()
                 if name.startswith('test_') and name != 'test_all' and
                 isinstance(obj, types.ModuleType)]
    for test in test_mods:
        test.main(show_plot=False)

if __name__ == '__main__':
    main()
