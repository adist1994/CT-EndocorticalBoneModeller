//
// Created by rap58 on 29/09/15.
//

#include <iostream>
#include "colourMap.h"

ColourMap::ColourMap() {
}

int ColourMap::getLength() {
    return length;
}

void ColourMap::getNanColour(double &r, double &g, double &b) {
    r= nanColour[0]; g= nanColour[1]; b= nanColour[2];
}

void ColourMap::getBelowColour(double &r, double &g, double &b) {
    r=below[0]; g=below[1]; b=below[2];
}

void ColourMap::getAboveColour(double &r, double &g, double &b) {
    r=above[0]; g=above[1]; b=above[2];
}

void ColourMap::getColour(int index, double &r, double &g, double &b) {
    r=map[index*3];
    g=map[index*3+1];
    b=map[index*3+2];
    if(index>length) {
        std::cerr<<"Error in ColourMap::getColour - length is too large "<<length<<std::endl;
    }
}


// define nan colour
double ColourMap::nanColour[3]= {150,150,150};
double ColourMap::below[3]=     {106,  2,  44};
double ColourMap::above[3]=     {67,  50, 106};

// generate the colourmap
double ColourMap::map[length*3]={
        159,   3,  67,
        161,   5,  68,
        163,   7,  69,
        165,   9,  70,
        167,  12,  72,
        169,  14,  73,
        172,  17,  74,
        174,  19,  75,
        176,  22,  76,
        178,  25,  77,
        181,  27,  78,
        183,  30,  78,
        185,  32,  79,
        188,  35,  80,
        190,  38,  80,
        193,  40,  81,
        195,  43,  81,
        197,  45,  81,
        199,  48,  81,
        202,  50,  81,
        204,  52,  81,
        206,  55,  81,
        208,  57,  80,
        210,  59,  80,
        212,  61,  79,
        214,  63,  78,
        216,  65,  78,
        218,  67,  77,
        220,  69,  76,
        222,  71,  75,
        223,  73,  74,
        225,  74,  73,
        226,  76,  72,
        228,  78,  71,
        229,  80,  71,
        231,  81,  70,
        232,  83,  69,
        233,  85,  68,
        234,  87,  68,
        235,  88,  67,
        236,  90,  66,
        237,  92,  66,
        238,  94,  66,
        239,  95,  66,
        240,  97,  65,
        241,  99,  65,
        242, 101,  66,
        242, 103,  66,
        243, 105,  66,
        244, 108,  67,
        244, 110,  67,
        245, 112,  68,
        245, 114,  69,
        246, 117,  70,
        246, 119,  71,
        247, 121,  72,
        247, 124,  73,
        248, 126,  74,
        248, 129,  75,
        248, 132,  77,
        249, 134,  78,
        249, 137,  79,
        250, 140,  81,
        250, 142,  82,
        250, 145,  83,
        251, 148,  85,
        251, 150,  86,
        251, 153,  88,
        251, 156,  89,
        252, 159,  90,
        252, 161,  92,
        252, 164,  93,
        252, 167,  94,
        253, 169,  95,
        253, 172,  97,
        253, 174,  98,
        253, 177,  99,
        253, 179, 100,
        253, 182, 101,
        253, 184, 102,
        253, 186, 103,
        254, 188, 104,
        254, 190, 105,
        254, 193, 107,
        254, 195, 108,
        254, 197, 109,
        254, 199, 111,
        254, 201, 112,
        254, 203, 114,
        254, 204, 115,
        254, 206, 117,
        254, 208, 119,
        254, 210, 121,
        254, 212, 123,
        254, 213, 125,
        254, 215, 127,
        254, 217, 130,
        254, 218, 132,
        254, 220, 135,
        254, 222, 138,
        254, 223, 140,
        254, 225, 143,
        254, 226, 146,
        254, 228, 149,
        254, 229, 152,
        254, 231, 155,
        254, 232, 158,
        254, 234, 161,
        254, 235, 164,
        254, 237, 166,
        254, 238, 169,
        254, 239, 172,
        254, 241, 174,
        254, 242, 177,
        254, 243, 179,
        255, 244, 181,
        255, 246, 183,
        255, 247, 185,
        255, 248, 187,
        255, 249, 188,
        255, 250, 189,
        255, 251, 190,
        255, 252, 190,
        255, 252, 191,
        255, 253, 191,
        255, 254, 191,
        255, 254, 191,
        254, 254, 190,
        254, 255, 189,
        254, 255, 188,
        253, 255, 187,
        253, 255, 186,
        252, 255, 184,
        252, 255, 182,
        251, 255, 181,
        250, 254, 179,
        249, 254, 177,
        248, 254, 175,
        247, 253, 172,
        246, 253, 170,
        245, 252, 168,
        244, 251, 166,
        243, 250, 164,
        241, 250, 162,
        240, 249, 160,
        238, 248, 159,
        237, 247, 157,
        235, 246, 155,
        233, 245, 154,
        231, 244, 153,
        229, 243, 152,
        227, 242, 151,
        225, 241, 150,
        223, 240, 149,
        221, 239, 149,
        219, 238, 149,
        216, 236, 149,
        214, 235, 149,
        212, 234, 149,
        209, 233, 149,
        207, 232, 150,
        204, 231, 151,
        202, 230, 151,
        200, 229, 152,
        197, 228, 153,
        195, 227, 154,
        192, 227, 155,
        190, 226, 156,
        187, 225, 157,
        185, 224, 158,
        183, 223, 159,
        180, 222, 161,
        178, 222, 162,
        175, 221, 162,
        173, 220, 163,
        170, 219, 164,
        168, 219, 165,
        165, 218, 165,
        163, 217, 166,
        160, 216, 166,
        158, 215, 167,
        155, 215, 167,
        153, 214, 167,
        150, 213, 167,
        147, 212, 167,
        145, 211, 167,
        142, 210, 167,
        139, 209, 167,
        137, 208, 167,
        134, 207, 166,
        131, 206, 166,
        128, 205, 166,
        125, 204, 166,
        123, 203, 165,
        120, 201, 165,
        117, 200, 165,
        114, 199, 165,
        111, 197, 165,
        108, 196, 165,
        105, 194, 165,
        102, 193, 165,
        99, 191, 165,
        95, 189, 166,
        92, 187, 166,
        89, 185, 167,
        86, 183, 167,
        83, 181, 168,
        80, 179, 169,
        78, 177, 170,
        75, 175, 171,
        72, 173, 172,
        70, 171, 173,
        67, 168, 174,
        65, 166, 176,
        63, 164, 177,
        61, 161, 178,
        59, 159, 179,
        57, 156, 181,
        55, 154, 182,
        54, 151, 183,
        53, 148, 184,
        52, 146, 186,
        51, 143, 187,
        51, 140, 187,
        50, 138, 188,
        50, 135, 189,
        50, 132, 190,
        50, 130, 190,
        51, 127, 190,
        51, 124, 190,
        52, 122, 190,
        53, 119, 190,
        54, 116, 190,
        56, 114, 189,
        57, 111, 188,
        59, 109, 187,
        61, 106, 186,
        63, 104, 185,
        65, 101, 184,
        67,  99, 182,
        69,  97, 181,
        72,  95, 179,
        74,  93, 177,
        77,  91, 176,
        79,  89, 174,
        82,  87, 172,
        84,  85, 170,
        87,  83, 168,
        89,  82, 166,
        91,  80, 164,
        94,  79, 162,
        96,  77, 160,
        98,  76, 159
};
