#include <iostream>
#include <vector>

class coord{

public:
    void set_coords(int cx, int cy);
    int x;
    int y;
};

void coord::set_coords(int cx, int cy)
{
    x = cx;
    y = cy;
}

int main()
{

    std::vector<coord> test;

    coord temp;

    temp.set_coords(1,2);
    test.push_back(temp);

    temp.set_coords(3,4);
    test.push_back(temp);

    test[0].x += 5;
    test[0].y -= 3;

    for (int i = 0; i < test.size(); i++)
        std::cout<<test[i].x<<"\t"<<test[i].y<<std::endl;



}