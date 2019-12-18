#include <iostream>
#include <vector>
#include "aux.h"

using namespace std;

void gprint(vector<int> &v) {
    cout << "[ ";
    for(int i = 0; i < v.size(); i++) {
        cout << v[i] << " ";
    }
    cout << "]" << endl;
}
