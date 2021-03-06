Questa cartella contiene gli **esercizi del NSL (Numerical Simulation Laboratory)** - A.A. 2019/2020 - *Simone Codispoti* - matricola 916834.

Alcune note sul lavoro svolto:

- I **Jupyter Notebook** sono scritti in inglese (con più che probabili errori!) e suddivisi in opportuni paragrafi per indicare la soluzione di una particolare istanza dell'esercizio. Ogni jupyter è nominato "ExerciseNum", con "Num" il numero dell'esercitazione;

- Ogni esercitazione contiene un singolo **makefile**, che compila "exerciseNum.cpp", il main code della esercitazione, e genera l'eseguibile "exerciseNum.exe"; le esercitazioni 4 - 6 - 7 sfruttano invece il codice già fornito e opportunamente modificato, quindi sono associate a eseguibili di nome diverso. Ogni makefile contiene un metodo "clean" che rimuove i file oggetto, gli output files e il seed.out;

- Gli **output di ogni esercitazione** sono stati organizzati in opportune sottocartelle, in modo da rendere chiara la struttura del codice (librerie e classi utilizzate). Per la maggior parte delle esercitazioni, compilando il codice ed eseguendo il file.exe si ottengono gli stessi output di quelli presenti nelle sottocartelle citate. Questo *non* accade per le esercitazioni 4 - 6 - 7, nelle quali si è dovuto cambiare frequentemente i parametri di input ed eseguire più volte il codice: le equilibrazioni sono state svolte "a mano", con una procedura descritta dettagliatamente nei rispettivi jupyter;

- Nelle esercitazioni 8 - 9 - 10 si è scelto di lavorare marcatamente *Object Oriented*, focalizzando l'attezione sulla creazione delle classi **VQMC.h/.cpp** e **TSP.h/.cpp**: il main code di queste esrcitazioni contiene i dati del problema, le configurazioni iniziali e l'algoritmo di ottimizzazione adottato, che possono essere modificati ogni volta, lasciando invece invariate le classi;

- Infine, oltre alle librerie specifiche per ogni esercitazione e alla libreria **"random.h"**, si sono sfruttate ampiamente le librerie **"posizione.h"**, **"funzione.h"** (l'implementazione delle quali nei rispettivi file .cpp è abbastanza autoesplicativa), e **"utilities.h"**, contenente funzioni template per la stampa di output e, soprattutto, il metodo `template <class T> void MC_Mean_Error(const vector <T>& input, vector <T>& average, vector <T>& error, const int n_step, const int n_cell)`, che implementa il **Data-Blocking** sul vettore di dati `input` con `n_step` MC steps totali e `n_cell` blocchi, salvando medie e errori progressivi nei vettori `average` e `error`.
