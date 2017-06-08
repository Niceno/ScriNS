"""
Run all test
"""

import timeit
import types

from scrins import tests

def main():
    """
    Docstring.
    """

    test_mods = [(name, obj) for name, obj in vars(tests).items()
                 if name.startswith('test_') and name != 'test_all' and
                 isinstance(obj, types.ModuleType)]
    total_time = 0
    for name, test in test_mods:

        print('#' * 40)
        print(name)
        print('#' * 40)
        start = timeit.default_timer()
        test.main(show_plot=False)
        duration = timeit.default_timer() - start
        total_time += duration
        print('run time {}: {:7.3f} s'.format(name, duration))
        print('total time: {:7.3f} s'.format(total_time))

if __name__ == '__main__':
    main()
