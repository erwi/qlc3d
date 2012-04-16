#ifndef SHAPEFUNCTION2D_H
#define SHAPEFUNCTION2D_H

//
//  STACK ALLOCATED SURFACE SHAPE FUNCTIONS
//




#define NGPS4 6
namespace SurfaceShapes
{
    const double gps4[NGPS4][2]= {
        {0.8168476, 0.09157621},
        {0.09157621,0.8168476},
        {0.09157621,0.09157621},
        {0.1081030, 0.4459485},
        {0.4459485, 0.1081030},
        {0.4459485, 0.4459485}};
}
class ShapeSurf4thOrder{
public:
    const unsigned int ngps;
    double w[NGPS4];
    double sh1[NGPS4][3];
    double sh1r[NGPS4][3];
    double sh1s[NGPS4][3];
    double sh1t[NGPS4][3];

    ShapeSurf4thOrder():
        ngps(NGPS4)

    {
    w[0] = 0.05497587;
    w[1] = 0.05497587;
    w[2] = 0.05497587;
    w[3] = 0.1116908;
    w[4] = 0.1116908;
    w[5] = 0.1116908;
        for(unsigned int i = 0 ; i < NGPS4 ; i++)
        {
            // P1 Shape functions
            sh1[i][0]=1-SurfaceShapes::gps4[i][0]-SurfaceShapes::gps4[i][1];
            sh1[i][1]=SurfaceShapes::gps4[i][0];
            sh1[i][2]=SurfaceShapes::gps4[i][1];

            // P1 Shape functions r-derivatives
            sh1r[i][0]=-1;
            sh1r[i][1]=1;
            sh1r[i][2]=0;

            // P1 Shape functions s-derivatives
            sh1s[i][0]=-1;
            sh1s[i][1]=0;
            sh1s[i][2]=1;

            //P1 Shape functions t-derivatives
            sh1t[i][0]=-1;
            sh1t[i][1]=0;
            sh1t[i][2]=0;

        }
    }
};


#endif // SHAPEFUNCTION2D_H
