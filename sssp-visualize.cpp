#include "opencv2/opencv.hpp"
#include "sssp_common.h"

#define cvCellSizeX 65
#define cvCellSizeY 15
#define cvFontSize 0.3

extern void showSSSP(float *matrixArray, pathElement *intermediateCostsAndPath)
{
    cv::Mat sspImage(MATRIX_SIZE * cvCellSizeX, MATRIX_SIZE * cvCellSizeY, CV_8UC3, cv::Scalar(0, 0, 0));

    for (int i = 0; i < MATRIX_SIZE; i++)
    {
        for (int j = 0; j < MATRIX_SIZE; j++)
        {
            cv::rectangle(sspImage, cv::Point(j * cvCellSizeX, i * cvCellSizeY), cv::Point((j + 1) * cvCellSizeX, (i + 1) * cvCellSizeY), cv::Scalar(255, 255, 255), 1);
            std::string val = std::to_string(matrixArray[i * MATRIX_SIZE + j]);
            cv::putText(sspImage, std::to_string(j) /*std::to_string((MATRIX_SIZE - 10)*(i * MATRIX_SIZE + j))*/, cv::Point(j * cvCellSizeX + 5, i * cvCellSizeY + cvCellSizeY/2 + cvCellSizeY/3), cv::FONT_HERSHEY_SIMPLEX, cvFontSize, cv::Scalar(0, 255, 255), 1);
        }

    }
    cv::imshow("SSSP", sspImage);
    cv::waitKey(0);
    
    
}