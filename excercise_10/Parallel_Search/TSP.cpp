#include "TSP.h"

// --- Individual istances --- //

Individual :: Individual(int complexity){       // Individual initialization

    for(int i=0; i<complexity; i++)
        m_ind.push_back(City());

    m_fit = 0;      // fitness 0 for an "empty" Individual
}

void Individual :: Print_DNA(){         // Simple visualization of the gene informations

    cout << "-------------------------------------------------------" << endl << endl;
    for(int i=0; i<m_ind.size(); i++){
        cout << "Gene number " << i+1 << ": " << endl;
        cout << "City name:       " << m_ind[i].Get_name() << endl;
        cout << "City position:   " << "x = " << m_ind[i].Get_pos().Get_X() << "  " 
        << "y = " <<  m_ind[i].Get_pos().Get_Y() << "  " << "z = " <<  m_ind[i].Get_pos().Get_Z() << endl << endl; 
    }
    cout << "-------------------------------------------------------" << endl << endl;
}

void Individual :: Print_DNA(string filename){        // Getting city coordinates in an output file

    ofstream out;
    out.open(filename);
    if(out.fail()){
        cerr << "Could not open" << filename << "!" << endl;
        exit(1);
    }

    for(int i=0; i<m_ind.size(); i++){
        out << m_ind[i].Get_pos().Get_X() << setw(12) << m_ind[i].Get_pos().Get_Y() << endl;
    }
}

bool Individual :: DNA_Corrupted(){

    /*
        Chek for the DNA of a TSP individual: we must fulffill two bonds:
          - every gene must be initialized (i.e. a "City" must be created, at least with a precise position);
          - every gene must be different from each other (at least by the different posistions ---> see operator == overload);
          - from the previuos point applied to every gene follows that every gene is not repeted in the DNA;
    */
   
    bool corrupted = 0;
    int size = m_ind.size();

    int ctr = 0;
    for(int i=0; i<size; i++)
        if(m_ind[i].Get_pos().Get_X() == 0 && m_ind[i].Get_pos().Get_Y() == 0) ctr++;

    if(ctr > 1){
        corrupted = 1;
        return corrupted;
    }

    for(int i=0; i<size-1; i++ ){
        for(int j=i+1; j<size; j++){
            if(m_ind[i] == m_ind[j] && m_ind[i] == m_ind[j]){
                corrupted = 1;
                return corrupted;
            }
        }
    }

    return corrupted;
}

void Individual :: Eval_Fitness(){            // Fitness evaluation using a quadratic distance cost function

    int size = m_ind.size();
    double fitness = 0;

    for(int i=0; i<(size-1); i++)
        fitness += m_ind[i].Get_pos().Norm_Quad_R2(m_ind[i+1].Get_pos());

    fitness += m_ind[size-1].Get_pos().Norm_Quad_R2(m_ind[0].Get_pos());

    m_fit = fitness;
}

// --- GA operators --- //

// Every operator ensure the preservation of the first gene, to reduce degeneration in the evaluation of the path

void Individual :: Swap_Mutation(Random& rand){      // Simple interchange between two genes

    int size = m_ind.size();
    int last = size - 1;
    int pos = int(rand.Rannyu(1, size));

    if(pos == last) swap(m_ind[last], m_ind[1]);
    else swap(m_ind[pos], m_ind[pos + 1]);
}

void Individual :: Per_Mutation(Random& rand){      // A variance of Swap_Mutation: interchange of a sequence of genes

    int size = m_ind.size();
    int last = size - 1;
    int range = last/2;
    int pos = int(rand.Rannyu(2, range));

    for(int i=pos; i>=1; i--)
        swap(m_ind[i], m_ind[i + pos]);
}

void Individual :: Inversion_Mutation(Random& rand){    // Inversion of the order of genes

    int size = m_ind.size();
    int pos = int(rand.Rannyu(3, size));
    int step = pos/2;

    for(int i=0; i<step; i++)
        swap(m_ind[pos - i], m_ind[i + 1]);
}

void Individual :: Shift_Mutation(Random& rand){        // Uniform shift of a sequence of genes

    int size = m_ind.size();
    int last = size - 1;
    int range = last/2;
    int pos = int(rand.Rannyu(1, range));

    Individual mutant(size);            // We safely use another Individual as support
    mutant.Set_Gene(m_ind[0],0);
    int index = 0;
    for(int i=1; i<=pos; i++){
        index = last - pos + i;
        mutant.Set_Gene(m_ind[index], i);
    }
    for(int i=1; i<=(last-pos); i++)
        mutant.Set_Gene(m_ind[i], pos + i);

    for(int i=1; i<m_ind.size(); i++)
        m_ind[i] = mutant.Get_Gene(i);
}

// --- Population operations --- //

vector <Individual> Generation_0(Individual progenitor, int size, Random& rand){

    vector <Individual> gen_0;

    if(progenitor.DNA_Corrupted()){
        cerr << "Could not initialize the first generation because the progenitor has corrupted DNA!" << endl;
        cerr << "DNA sequence : " << endl << endl;
        progenitor.Print_DNA();
        cerr << "Exit program" << endl;
        exit(-1);
    }

    for(int i=0; i<size; i++){
        for(int j=0; j<3; j++)
            progenitor.Swap_Mutation(rand);       // we generate random individuals simply applying 3 Swap Mutations to the "progenitor"
    gen_0.push_back(progenitor);
    }

    return gen_0;
}

void Pop_Sorting(vector <Individual>& pop){          // population fitness evaluation adn sorting

    for(int i=0; i<pop.size(); i++)     // fitness evaluation of each individual
        pop[i].Eval_Fitness();

    sort(pop.begin(), pop.end());       // sorting from most fit to least fit ---> see operator < overload
}

Individual Natural_Selection(const vector <Individual>& pop, Random& rand){        // select an individual from the SORTED population with a probability proportional to its fitness

    int bias = 4;           // increase this parameter to obtain a stronger selection
    int N = pop.size();
    int pos = int(N*pow(rand.Rannyu(),bias));
    return pop[pos];
}

void Crossover(Individual& sel_1, Individual& sel_2, Random& rand){

    int size = sel_1.Get_Complexity();
    int last = size - 1;
    int pos = int(rand.Rannyu(1, size));

    int ctr = 0;
    Individual first_child;             // inheritance of the genes of the sel_2 ind different from the tail of the sel_1 ind
    first_child.Add_Gene(sel_2.Get_Gene(0));
    for(int i=1; i<=last; i++){
        ctr = 0;
        for(int j=pos; j<=last; j++){
            if(sel_2.Get_Gene(i) == sel_1.Get_Gene(j))
                ctr ++;
        }
        if(ctr==0) first_child.Add_Gene(sel_2.Get_Gene(i));
    }
    for(int l=pos; l<=last; l++)        // inheritance of the tail of the sel_1 ind
        first_child.Add_Gene(sel_1.Get_Gene(l));

    Individual second_child;             // inheritance of the genes of the sel_1 ind different from the tail of the sel_2 ind
    second_child.Add_Gene(sel_1.Get_Gene(0));
    for(int i=1; i<=last; i++){
        ctr = 0;
        for(int j=pos; j<=last; j++){
            if(sel_1.Get_Gene(i) == sel_2.Get_Gene(j))
             ctr ++;
        }
        if(ctr==0) second_child.Add_Gene(sel_1.Get_Gene(i));
    }
    for(int l=pos; l<=last; l++)        // inheritance of the tail of the sel_2 ind
        second_child.Add_Gene(sel_2.Get_Gene(l));

    sel_1 = first_child;
    sel_2 = second_child;

    if(sel_1.DNA_Corrupted() or sel_2.DNA_Corrupted()){
        cerr << "Crossover mutation failed!" << endl;
        cerr << "DNA sequences corrupted sequences: " << endl << endl;
        sel_1.Print_DNA();
        sel_2.Print_DNA();
        cerr << "Exit program" << endl;
        exit(-2);
    }
}