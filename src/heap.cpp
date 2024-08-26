#include "../lib/lib.h"

int getLeft(int parent)
{
    return parent * 2 + 1;
}

int getRight(int parent)
{
    return parent * 2 + 2;
}

int getParent(int child)
{
    return (child - 1) / 2;
}

void FMM::DeleteHeap(vector<vector<int>>& H, int &size)
{
    int left; 
    int right;
    int smaller; 
    int parent; 
    
    swap(H.at(0), H.at(size-1));

    size--;
    parent = 0;
    
    while(1){
        left = getLeft(parent);
        right = getRight(parent);
 
        if(left < size && right < size){
            if (H.at(left).at(0) > H.at(right).at(0)) {
                smaller = right;
            } else {
                smaller = left;
            }
        }else if(left < size){
            smaller = left;
        }else {
            break;
        }

        if(H.at(smaller).at(0) >= H.at(parent).at(0)){
            break;
        }
        swap(H.at(smaller), H.at(parent));

        parent = smaller;
    }
}

void FMM::UpHeap(vector<vector<int>> &H, int _i, int _j)
{
    int parent;
    int target;

    auto it = find_if(
        begin(H), end(H),
        [&_i, &_j](const auto& row) {
            return row.at(1) == _i && row.at(2) == _j;
        }
    );
    if(it == end(H)){
        cout << "not found";
    }else{  
        target = distance(begin(H), it);
    }

    while(1){    
        parent = getParent(target);

        if(H.at(parent).at(0) > H.at(target).at(0)){
            swap(H.at(parent), H.at(target));
            target = parent; 

            if (target == 0){
                break;
            }
        }else {
            break;
        }
    }
    
}