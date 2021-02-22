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

// Finds ensemble average of trace mobility matricies of i^th j^th beads
// that in 'ensembleTraceMobility[i][j]'. For why we need this see:
// "Diffusion coefficients of elastic macromolecules Bogdan Cichocki, 
// Marcin Rubin, Anna Niedzwiecka, and Piotr Szymczak; J. Fluid Mech. (2019)."
//
// Uses Rotne-Prager approximation which simplifies things considerably. 
// For mathematical expressions compare with: "Rotne-Prager-Yamakawa 
// approximation for different-sized particles in application to macromolecular bead
// models Pawel Zuk, Eligiusz Wajnryb, Krzysztof Mizerski, and Piotr Szymczak; 
// J. Fluid Mech. (2014)."

#pragma once
#include "chain_generator.cpp"
#include<map>
#include<cmath>

typedef std::map<int,std::map<int,double> > vararr;
typedef std::vector<Chain> Ensemble;

vararr ensembleAveragedMobility(Ensemble ensemble, vdouble ensembleBeadSizesDesign)
{
    int l = ensembleBeadSizesDesign.size();
    int m = ensemble.size();
    auto ensembleTraceMobility = vararr();
    for(int i=0;i<l;i++)for(int j=0;j<l;j++) ensembleTraceMobility[i][j] = 0;
    for(int i=0;i<l;i++)
    {
        for(int j=0;j<l;j++)
        {
            if(i==j)
            {
                // Self mobility
                ensembleTraceMobility[i][j] = (3.0/(6.0*M_PI*ensembleBeadSizesDesign[i]));
            }
            else
            {
                double ai = ensembleBeadSizesDesign[i];
                double aj = ensembleBeadSizesDesign[j];
                double accum = 0.0;
                // enseble average
                for(int k=0;k<m;k++)
                {
                    double Rij = distance( ensemble[k][i] , ensemble[k][j] );
                    if(Rij > ai + aj) // far appart
                    {
                        accum = accum + (3.0/(6.0*M_PI*Rij));
                    }
                    else // overlapping
                    {
                        double TTidentityScale = (1.0 / (6.0 * M_PI * ai * aj)) * ((16.0*(pow(Rij,3))*(ai+aj)-pow(pow(ai-aj,2) + 3*(pow(Rij,2)),2))/(32.0 * (pow(Rij,3))));
                        double TTrHatScale = (1.0 / (6.0 * M_PI * ai * aj)) * ( 3 * pow(pow((ai-aj),2) - pow(Rij,2),2) / (32 * (pow(Rij,3))) );

                        accum = accum + TTrHatScale + 3*TTidentityScale;
                    }
                }
                accum = accum / m;
                ensembleTraceMobility[i][j] = accum;
            }
        }
    }
    return ensembleTraceMobility;
}
