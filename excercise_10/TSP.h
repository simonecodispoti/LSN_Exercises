#ifndef TSP_H
#define TSP_H
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include "posizione.h"
#include "random.h"
using namespace std;

//  || --- | City class | --- ||  //

class City{

    public:

        // Constructors
        City(){}
        City(Posizione pos){m_pos = pos;}
        City(Posizione pos, string name){m_pos = pos;
                                        m_name = name;}
        ~City(){};

        // Methods
        void Set_name(string name) {m_name = name;}
        string Get_name() const {return m_name;}
        void Set_pos(Posizione pos) {m_pos = pos;}
        Posizione Get_pos() const {return m_pos;}

        // Operators
        bool operator ==(const City& rhs){      // Two cities are identical if they are in the same position 
            if(m_pos.Get_X() == rhs.Get_pos().Get_X() && m_pos.Get_Y() == rhs.Get_pos().Get_Y())
                return 1;
            else
                return 0;
        }

    private:

        string m_name;
        Posizione m_pos;
};


//  || --- | TSP class individuals with GA operators | --- ||  //

class Individual{

    public:

        // Constructors
        Individual(){}
        Individual(int complexity);
        ~Individual(){}

        // Methods
        int Get_Complexity() const {return m_ind.size();}
        City Get_Gene(int position) const {return m_ind[position];}
        void Set_Gene(City city, int position) {m_ind[position] = city;}
        void Add_Gene(City city) {m_ind.push_back(city);}

        // Operators
        bool operator <(const Individual& rhs) {return m_fit < rhs.Get_Fitness();}

        // Gene visualization and control
        void Print_DNA();
        void Print_DNA(string filename);
        bool DNA_Corrupted();

        // Fitness function
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

// Crossover 
void Crossover(Individual& sel_1, Individual& sel_2, Random& rand);

//  || --- | Population operations | --- ||  //

vector <Individual> Generation_0(Individual progenitor, int size, Random& rand);
void Pop_Sorting(vector <Individual>& pop);
int Natural_Selection(const vector <Individual>& pop, Random& rand);

#endif /* TSP_H */