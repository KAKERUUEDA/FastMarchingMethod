#include "lib.h"

int getLeft(int parent) {
    return parent * 2 + 1;
}

int getRight(int parent) {
    return parent * 2 + 2;
}

int getParent(int child) {
    return (child - 1) / 2;
}
/*
void FMM::AddHeap(vector<int>& H, int size) {

    int add; 
    int parent;

    add = size;
    if (add == 0) {
        return;
    }

    while (1) {
        parent = getParent(add);
        
        if (H.at(parent) > H.at(add)) {
            swap(H.at(parent), H.at(add));
            add = parent;
            if (add == 0) {
                break;
            }
        } else {
            break;
        }
    }
}*/
void FMM::DeleteHeap(vector<vector<int>>& H, int &size) {
    int left; 
    int right;
    int smaller; 
    int parent; 
       //cout << "前の葉のi座標は = " << H.at(size-1).at(1) <<  " 葉のj座標は = " << H.at(size-1).at(2) << endl;
    swap(H.at(0), H.at(size-1));

    //cout << "size = " << size << endl;
    //cout << "新しい根のi座標は = " << H.at(0).at(1) << " 根のj座標は = " << H.at(0).at(2) <<  endl;
    //cout << "消された葉のi座標は = " << H.at(size-1).at(1) <<  " 葉のj座標は = " << H.at(size-1).at(2) << endl;
    
    
    //H_tmp.at(tmp2).push_back(H.at(size-1).at(0));
    //H_tmp.at(tmp2).push_back(H.at(size-1).at(1));
    //H_tmp.at(tmp2).push_back(H.at(size-1).at(2));

    //tmp2++;

    size--;
    //cout << size << endl;
    parent = 0;
    
    while (1) {
        left = getLeft(parent);
        right = getRight(parent);
        //cout << "left = " << left  << " right = " << right << endl;
 
        if (left < size && right < size) {
            if (H.at(left).at(0) > H.at(right).at(0)) {
                smaller = right;
            } else {
                smaller = left;
            }
        } else if (left < size) {
            smaller = left;
        } else {
            break;
        }

        if (H.at(smaller).at(0) >= H.at(parent).at(0)) {
            break;
        }
        swap(H.at(smaller), H.at(parent));

        parent = smaller;
        //cout << "parent = " << parent << endl;
    }
}
/*
void makeHeap(int a[], int num) {

    int size;
    size = 0;

    while (size < num) {

        addHeap(a, size);
        size++;
    }
}*/



void FMM::UpHeap(vector<vector<int>> &H, int _i, int _j){
    int parent;
    int target;


    auto it = find_if(
        begin(H), end(H),
        [&_i, &_j](const auto& row) {
            return row.at(1) == _i && row.at(2) == _j;
        }
    );
    //cout << "2_1_6_1" << endl;
    if(it == end(H)) {
        cout << "not found";
    }
    else {  
        //cout << "nearの位置は根の0番目から数えて" << distance(begin(H), it) << endl;
        target = distance(begin(H), it);
        }

    while (1) {    

        parent = getParent(target);
        //cout << "parent = " << parent << endl;
        if (H.at(parent).at(0) > H.at(target).at(0)) {

            swap(H.at(parent), H.at(target));
            target = parent; 

            if (target == 0) {
                break;
            }
                   //cout << "2_1_6_3" << endl;
        }else {
            break;
        }
    }
    
}