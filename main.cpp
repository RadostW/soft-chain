//Copyright 2021 Radost Waszkiewicz
//Permission is hereby granted, free of charge, to any person obtaining a copy 
//of this software and associated documentation files (the "Software"), to deal 
//in the Software without restriction, including without limitation the rights 
//to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
//copies of the Software, and to permit persons to whom the Software is 
//furnished to do so, subject to the following conditions:
//The above copyright notice and this permission notice shall be included in all
//copies or substantial portions of the Software.
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
//IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
//FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
//AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
//LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE 
//SOFTWARE.

#include "chain_generator.cpp"
#include "rotne_prager.cpp"
#include "args_parser.cpp"
#include "Eigen/Dense"
#include <stdio.h>
#include <iostream>

int main(int argc, char **argv)
{
    InputParser input(argc, argv);
    if(input.cmdOptionExists("-h") || input.cmdOptionExists("--help")){
        printf("Soft chain generator!\n\n");
        printf("Chain descriptions are json arrays of triples [hydrodynamic size, steric size, count]\n");
        printf("\n");
        printf("Usage:\n");
        printf("main.exe [-h|--help]\n");
        printf("main.exe -c CHAIN -n ENSEMBLESIZE [-m|-r]\n");
        printf("\n");
        printf("Sample single chain of length 10:\nmain.exe -c [[1,1,10]]\n\n");
        printf("Sample single chain with variable beadsizes:\nmain.exe -c [[1,1,1],[2,2,1],[3,3,1]]\n\n");
        printf("Sample 33 chains of length 10:\nmain.exe -c [[1,1,10]] -n 33\n\n");
        printf("Grand mobility matrix from ensemble of size 33:\nmain.exe -c [[1,1,10]] -n 33 -m\n\n");
        printf("Hydrodynamic size from ensemble of size 33:\nmain.exe -c [[1,1,10]] -n 33 -r\n\n");
    }
    else
    {
        if(input.cmdOptionExists("-c"))
        {
            auto buffer = input.getCmdOption("-c");
            std::vector<double> hydrodynamicSizes;
            std::vector<double> stericSizes;
            double h,s;
            int n;
            std::replace( buffer.begin(), buffer.end(), '[', ' ');
            std::replace( buffer.begin(), buffer.end(), ']', ' ');
            std::replace( buffer.begin(), buffer.end(), ',', ' ');
            int consumed=0,chars=0;
            while(sscanf(buffer.c_str()+consumed,"%lf %lf %d %n",&h,&s,&n,&chars) == 3)
            {
                consumed +=chars; //Update number of chars consumed by sscanf.
                for(int i=0;i<n;i++)
                {
                    hydrodynamicSizes.push_back(h);
                    stericSizes.push_back(s);
                }
            }

            int ensembleSize = 1;
            if(input.cmdOptionExists("-n"))
            {
                sscanf(input.getCmdOption("-n").c_str(),"%d",&ensembleSize);
            }

            auto ensemble = Ensemble();
            for(int i=0;i<ensembleSize;i++)
            {
                ensemble.push_back(getChain(stericSizes.size(),0,stericSizes));
            }

            auto mobilityMatrix = ensembleAveragedMobility(ensemble,hydrodynamicSizes);

            if(input.cmdOptionExists("-m"))
            {
                //Output matrix
                printf("[");
                for(int i=0;i<mobilityMatrix.size();i++)
                {
                    printf("[");
                    for(int j=0;j<mobilityMatrix[0].size();j++)
                    {
                        printf("%lf",mobilityMatrix[i][j]);
                        if(j+1<mobilityMatrix[0].size())printf(",");
                    }
                    printf("]");
                    if(i+1<mobilityMatrix.size())printf(",");
                }
                printf("]");
            }
            else if(input.cmdOptionExists("-r"))
            {
                //Output hydrodynamic size
                int size = mobilityMatrix.size();
                Eigen::MatrixXd mobilityMatrixEigenized(size,size); 
                for(int i=0;i<size;i++)for(int j=0;j<size;j++)mobilityMatrixEigenized(i,j) = mobilityMatrix[i][j];

                std::cout <<"{";
                std::cout << "\"chain\": "<<input.getCmdOption("-c") << ", ";
                std::cout << "\"Cichocki\": " << mobilityMatrixEigenized.inverse().sum() / (2*M_PI) << ", ";
                std::cout << "\"Kirkwood\": " << (1.0/(2.0*M_PI)) * (1.0 * size * size / mobilityMatrixEigenized.sum()) << "}\n";
                
            }
            else
            {
                //Output sampled chains
                printf("[");
                for(int i=0;i<ensembleSize;i++)
                {
                    printf("[");
                    for(int j=0;j<stericSizes.size();j++)
                    {
                        printf("[[%lf,%lf,%lf],%lf]",ensemble[i][j].x,ensemble[i][j].y,ensemble[i][j].z,hydrodynamicSizes[j]);
                        if(j+1<stericSizes.size())
                        {
                            printf(",");
                        }
                    }
                    if(i+1<ensembleSize)
                    {
                        printf("],");
                    }
                    else
                    {
                        printf("]");
                    }
                }
                printf("]\n");
            }
        }
        else
        {
            printf("No chain description. Try --help.");
        }
    }
}
