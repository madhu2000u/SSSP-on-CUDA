#include "opencv2/opencv.hpp"
#include "sssp_common.h"

#define cvCellSizeX 65
#define cvCellSizeY 15
#define cvFontSize 0.6

float traversePath(cv::Mat image, int recursionCount, float *d_c, int index, pathElement *path)
{   
    int nextIndex;
    float cost;
    // int recursionCount;
    if(path[index].pathIndex == -1)
    {   
        printf("(%d, %d) -> ", (int)(index) / MATRIX_SIZE, index % MATRIX_SIZE);
        image.at<cv::Vec3b>((int)(index / MATRIX_SIZE)*0.9, (index % MATRIX_SIZE)*1) = cv::Vec3b(0, 255, 255);
        return d_c[index];
    };

    nextIndex = path[index].pathIndex;
    recursionCount++;
    float val = traversePath(image, recursionCount, d_c, nextIndex, path);
    recursionCount--;
    cost = (recursionCount > 0 ? d_c[index] : 0) + val;                                             //the first calling function is the one with destination whose weight should not be added to cost.
    printf("(%d, %d) -> ", (int)(index / MATRIX_SIZE), index % MATRIX_SIZE);
    image.at<cv::Vec3b>((int)(index / MATRIX_SIZE)*0.9, (index % MATRIX_SIZE)*1) = cv::Vec3b(0, 255, 255);
    // cv::circle(image, cv::Point(((nextIndex / MATRIX_SIZE)), ((int)(index % MATRIX_SIZE))), 1, cv::Scalar(0, 0, 255), -1);
    return cost;
}


extern void showSSSP(float *d_c, pathElement *intermediateCostsAndPath, int destinationIndex)
{
    // cv::Mat sspImage = cv::Mat::zeros(MATRIX_SIZE, MATRIX_SIZE, CV_8UC3);
     cv::Mat sspImage = cv::Mat::zeros(MATRIX_SIZE * 0.9, MATRIX_SIZE*1, CV_8UC3);
    

//     for (int i = 0; i < MATRIX_SIZE; i++)
//     {
//         for (int j = 0; j < MATRIX_SIZE; j++)
//         {
            
//             // cv::circle(sspImage, cv::Point(j * 1.5, i * 1.5), 0.5, cv::Scalar(0, 255, 255), -1);
//             // cv::rectangle(sspImage, cv::Point(j * cvCellSizeX, i * cvCellSizeY), cv::Point((j + 1) * cvCellSizeX, (i + 1) * cvCellSizeY), cv::Scalar(255, 255, 255), 1);
//             // std::string val = std::to_string(matrixArray[i * MATRIX_SIZE + j]);
//             // if(i*MATRIX_SIZE+j == 1023)cv::putText(sspImage, std::to_string(j) /*std::to_string((MATRIX_SIZE - 10)*(i * MATRIX_SIZE + j))*/, cv::Point(j * 4, i * 4), cv::FONT_HERSHEY_SIMPLEX, cvFontSize, cv::Scalar(0, 255, 255), 1);
//         }

//     }


    printf("\n");
    int recursionCount = 0;
    float cost = traversePath(sspImage, recursionCount, d_c, destinationIndex, intermediateCostsAndPath);
    printf("\ncalculated cost - %f", cost);
    cv::imshow("SSSP", sspImage);
    cv::waitKey(0);
return;

    // /// test
    // int width = 1024;
    // int height = 1024;

    // // Create a black image
    // cv::Mat image = cv::Mat::zeros(height * 0.8, width * 10, CV_8UC3);

    // // Calculate the number of points in each dimension
    // int numPointsX = width / 1;
    // int numPointsY = height / 1;

    // // Generate evenly spaced points for demonstration purposes
    // std::vector<cv::Point> points;
    // for (int y = 0; y < numPointsY; ++y) {
    //     for (int x = 0; x < numPointsX; ++x) {
    //         int pointX = x * 4;
    //         int pointY = y * 0.8;
    //         points.push_back(cv::Point(pointX, pointY));
    //     }
    // }

    // // Set the color of the points (in BGR format, e.g., CV_RGB(0, 0, 255) for red)
    // cv::Scalar color(0, 0, 255); // Red color

    // // Draw all points on the image
    // for (const cv::Point& point : points) {
    //     // Set the color of the pixel at the specified point
    //     image.at<cv::Vec3b>(point.y, point.x) = cv::Vec3b(color[0], color[1], color[2]);
    // }

    // // Display the image
    // cv::imshow("Points Example", image);
    // cv::waitKey(0);

    // // Check the number of points generated
    // printf("point size - %d", points.size());

    
    
}