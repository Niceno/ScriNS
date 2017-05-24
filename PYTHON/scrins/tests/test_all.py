"""
Run all test
"""

from tests import test_advection_1D, test_inlet_outlet


def main():
    tests = [test_advection_1D, test_inlet_outlet]
    for test in tests:
        test.main(show_plot=False)

if __name__ == '__main__':
    main()
