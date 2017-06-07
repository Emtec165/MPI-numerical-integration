#include "mpi.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>

const int master = 0;
const int intervalBeg = 0;
const int intervalEnd = 1;


double evaluateFOfX(int degree, double coeffs[], double x){
    double sum = 0.0;
    double powOfX;
    for (int i = 0; i < degree; i++) {
        powOfX = pow(x, i);//in first loop i == 0, then x = 1, then coeffs[0] * 1 == free expression
        sum += powOfX * coeffs[i];
    }
    return sum;
}


double analyticalFOfX(int degree, double coeffs[], double x){
    double sum = 0.0;
    double singleFraction;
    double powOfX;
    for (int i = 0; i < degree; i++) {
        powOfX = pow(x, i + 1);//factor + 1 (x^1 instead of x^0)
        singleFraction = powOfX * coeffs[i];
        singleFraction /= i + 1;//factor + 1 (x/1 instead od x/0)
        sum += singleFraction;
    }
    return sum;
}

int main(int argc, char *argv[]) {
    int rank, size;
    std::string fileName;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);


    int degree = 0, integration = 0;
    double *coeffs = NULL, *interval = new double[2];

    char *packBuffer = NULL;
    int packSize = 0, packPosition = 0;




    if (argc == 2){
        fileName = argv[1];
//        fileName = "test4.in";
    } else {
        std::cout << "Could not read file name. Aborting\n";
        MPI_Finalize();
        return 1;
    }

    if (rank == master){
        std::string readLine;
        std::ifstream inputFile;
        inputFile.open(fileName.c_str());


        if (inputFile.good()){

            int tmpPackSize;

            std::stringstream ss;
            getline(inputFile, readLine);
            degree = std::stoi(readLine.substr(6));//to skip "degree" letters
            degree += 1;//+ 1 because there is 0 factor degree

            coeffs = new double[degree];
            getline(inputFile, readLine);
            ss.str(readLine.substr(6));//6 to skip "coeffs" letters
            for (int i = 0; i < degree; i++) {
                ss >> coeffs[i];
            }

            getline(inputFile, readLine);
            ss = std::stringstream();//re-declare ss to clean ss
            ss.str(readLine.substr(8));//8 to skip "interval" letters
            ss >> interval[intervalBeg];
            ss >> interval[intervalEnd];

            getline(inputFile, readLine);
            integration = std::stoi(readLine.substr(11));//11 to skip "integration" letters





            MPI_Pack_size(2, MPI_INT, MPI_COMM_WORLD, &tmpPackSize);//[one int] coefficients count + [one int] integration (how many sectors)
            packSize = tmpPackSize;
            MPI_Pack_size(degree + 2, MPI_DOUBLE, MPI_COMM_WORLD, &tmpPackSize);//coefficients + boundary beg & end
            packSize += tmpPackSize;

            packBuffer = new char[packSize];

            MPI_Pack(&degree, 1, MPI_INT, packBuffer, packSize, &packPosition, MPI_COMM_WORLD);
            MPI_Pack(&coeffs[0], degree, MPI_DOUBLE, packBuffer, packSize, &packPosition, MPI_COMM_WORLD);
            MPI_Pack(&interval[0], 2, MPI_DOUBLE, packBuffer, packSize, &packPosition, MPI_COMM_WORLD);
            MPI_Pack(&integration, 1, MPI_INT, packBuffer, packSize, &packPosition, MPI_COMM_WORLD);


            MPI_Bcast(&packSize, 1, MPI_INT, master, MPI_COMM_WORLD);
            MPI_Bcast(packBuffer, packSize, MPI_PACKED, master, MPI_COMM_WORLD);

        } else {
            std::cout << "Unable to open file. Aborting\n";
            MPI_Finalize();
            return 1;
        }
        inputFile.close();

    } else {// slave

        MPI_Bcast(&packSize, 1, MPI_INT, master, MPI_COMM_WORLD);
        packBuffer = new char[packSize];
        MPI_Bcast(packBuffer, packSize, MPI_PACKED, master, MPI_COMM_WORLD);

        MPI_Unpack(packBuffer, packSize, &packPosition, &degree, 1, MPI_INT, MPI_COMM_WORLD);
        coeffs = new double[degree];
        MPI_Unpack(packBuffer, packSize, &packPosition, &coeffs[0], degree, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack(packBuffer, packSize, &packPosition, &interval[0], 2, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Unpack(packBuffer, packSize, &packPosition, &integration, 1, MPI_INT, MPI_COMM_WORLD);
    }

    std::cout << "Rank: " << rank;
    std::cout << " packSize: " << packSize;
    std::cout << " degree: " << degree;
    std::cout << " coeffs: ";
    for (int i = 0; i < degree; i++) {
        std::cout << coeffs[i] << " ";
    }
    std::cout << " interval: " << interval[intervalBeg] << " " << interval[intervalEnd];
    std::cout << " integration: " << integration << std::endl;






    int p = integration / size;
    int r = integration % size;
    int numberOfSubintervals, firstMidpoint;

    if (rank < r){
        numberOfSubintervals = p + 1;
        firstMidpoint = rank * (p + 1);
    } else {
        numberOfSubintervals = p;
        firstMidpoint = r * (p + 1) + (rank - r) * p;
    }



    double height = (interval[intervalEnd] - interval[intervalBeg]) / integration;
    double lOfR = interval[intervalBeg] + firstMidpoint * height;
    double uOfR = lOfR + numberOfSubintervals * height;
    std::cout << "Rank: " << rank << " l(r): " << lOfR << " u(r): " << uOfR << std::endl;





    double sum = 0.0;
    double numericalEvaluation;
    double analyticalEvaluation;
    for (int i = 0; i < numberOfSubintervals; i++) {
        double x1 = lOfR + i * height;
        double x2 = lOfR + (i + 1) * height;
        sum += evaluateFOfX(degree, coeffs, x1) + evaluateFOfX(degree, coeffs, x2);
    }
    sum /= 2;
    sum *= height;


    MPI_Reduce(&sum, &numericalEvaluation, 1, MPI_DOUBLE, MPI_SUM, master, MPI_COMM_WORLD);
    if (rank == master) {
        std::cout << "Numerical: " << numericalEvaluation << std::endl;

        analyticalEvaluation = analyticalFOfX(degree, coeffs, interval[intervalEnd]) - analyticalFOfX(degree, coeffs, interval[intervalBeg]);
        std::cout << "Analytical: " << analyticalEvaluation << std::endl;
    }






    MPI_Finalize();
    return 0;
}