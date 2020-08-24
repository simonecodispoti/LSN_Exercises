#ifndef TSP_H
#define TSP_H
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "posizione.h"
#include "random.h"
using namespace std;

//  || --- | City class | --- ||  //

class City{

    public:

        City(){}
        City(Posizione pos){m_pos = pos;}
        City(Posizione pos, string name){m_pos = pos;
                                        m_name = name;}
        ~City(){};

        void Set_name(string name) {m_name = name;}
        string Get_name() const {return m_name;}
        void Set_pos(Posizione pos) {m_pos = pos;}
        Posizione Get_pos() const {return m_pos;}

    private:

        string m_name;
        Posizione m_pos;
};


//  || --- | TSP class individuals with GA operators | --- ||  //

class Individual{

    public:

        Individual(){}
        Individual(int complexity);
        ~Individual(){}

        int Get_Complexity() const {return m_ind.size();}
        City Get_Gene(int position) const {return m_ind[position];}
        void Set_Gene(City city, int position) {m_ind[position] = city;}

        void Print_DNA();
        void Print_DNA(string filename);
        bool DNA_Corrupted();

        double Get_Fitness() const {return m_fit;};
        void Eval_Fitness();

        //  || --- | GAO (Generic Algorithm operators) | --- ||  //

        /* 
            These operators are implemented for this specific problem, but can be easily translated into other problems
            changing the individual type: instead of a "City" the data structure can be anything: "double", "binary", ...
        */

        void Swap_Mutation(Random& rand);
        void Per_Mutation(Random& rand);        
        void Inversion_Mutation(Random& rand);
        void Shift_Mutation(Random& rand);

    private:

        vector <City> m_ind;
        double m_fit;
};

//  || --- | Population operations | --- ||  //

vector <Individual> Generation_0(Individual progenitor, int size, Random& rand);
int Natural_Selection(const vector <Individual>& pop, Random& rand);
void Crossover(vector <Individual>& pop, int sel_1, int sel_2, Random& rand);

#endif /* TSP_H */