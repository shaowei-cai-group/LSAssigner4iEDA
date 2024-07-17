#include "test.h"
#include <iostream>
using namespace std;
int main(int argc, char *argv[])
{
    cout << "before test_panel" << endl;
    test_panel(argc, argv);
    cout << "after test_panel" << endl;
    return 0;
}