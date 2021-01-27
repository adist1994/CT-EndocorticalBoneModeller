//
// Created by rap58 on 29/09/15.
//

#ifndef MYPROJECT_COLOURMAP_H
#define MYPROJECT_COLOURMAP_H


class ColourMap {

public:
    ColourMap(void);

    int getLength(); // number of colours
    void getNanColour(double &r, double &g, double &b);
    void getBelowColour(double &r, double &g, double &b);
    void getAboveColour(double &r, double &g, double &b);
    void getColour(int index, double &r, double &g, double &b);


private:
    static constexpr int length = 253;

    static double map[3*length];
    static double nanColour[3];
    static double below[3];
    static double above[3];

};


#endif //MYPROJECT_COLOURMAP_H
