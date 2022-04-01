// Repo github: https://github.com/deeaanghelache/Genetic-Algorithms

/*
    Andreea Anghelache
    Grupa 251
    Tema Algoritmi Genetici
*/

// Genetic Algorithms

#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <random>
#include <stack>

using namespace std;

int dimensiuneaPopulatiei, capatStangaInterval, capatDreaptaInterval, a, b, c, precizie, numarEtape;
int lungimeCromozom;
double probabilitateRecombinare, probabilitateMutatie;

class Individ {
    double fitnessIndivid;
    double valoareIndivid;
    int valoareIntreaga;
    std::vector<int> cromozom;
public:
//    Individ(double valoareIndivid);
//    Individ(double valoareIndivid, int valoareIntreaga);
    Individ(double fitnessIndivid, double valoareIndivid, int valoareIntreaga, const vector<int> &cromozom);

    // Getters si setters
    double getFitnessIndivid() const;
    void setFitnessIndivid(double fitnessIndivid);
    double getValoareIndivid() const;
    void setValoareIndivid(double valoareIndivid);
    const std::vector<int> &getCromozom() const;
    void setCromozom(const std::vector<int> &cromozom);
    int getValoareIntreaga() const;
    void setValoareIntreaga(int valoareIntreaga);
};

/* Functii ajutatoare pentru algoritm */

double fitness(double x, const int a1 = -1, const int b1 = 1, const int c1 = 2){
    return a1 * x * x + b1 * x + c1;
}

vector<int> conversieZecimalInBinar(int valoare){
    vector<int> cromozom;
    stack<int> myStack;

    while(valoare > 0){
        myStack.push(valoare % 2);
        valoare /= 2;
    }

    int diferenta = lungimeCromozom - myStack.size();

    // ma asigur ca am un cromozom de size = lungimeCromozom (daca e nevoie, adaug 0-uri in fata)
    while(diferenta > 0){
        myStack.push(0);
        diferenta--;
    }

    while(!myStack.empty()){
        cromozom.push_back(myStack.top());
        myStack.pop();
    }

    return cromozom;
}

int putere10(int precizieInt){
    int putere = 1;
    for(int i = 0; i < precizieInt; i++){
        putere *= 10;
    }
    return putere;
}

int calculLungimeCromozom(){
    // log2((b-a) * 10^p)
    int putere = putere10(precizie);
    int lgCromozom = ceil(log2((b - a) * putere));

    return lgCromozom;
}

int genereazaValoareIntreaga(){
    // Translatarea intervalului

    int copieCapatStanga = capatStangaInterval, copieCapatDreapta = capatDreaptaInterval, putere = putere10(precizie);
    copieCapatStanga *= putere; // capat * 10^p
    copieCapatDreapta *= putere;

    copieCapatDreapta -= copieCapatStanga; // in exemplu: 2 * 10^6 - (-1) * 10^6
    copieCapatStanga -= copieCapatStanga; // 0

    // poate sa genereze numere care se scriu pe mai mult de lungimeCromozom biti
    int capatDreapta = min(copieCapatDreapta, (int) pow(2, lungimeCromozom) - 1);

    // Generare numere din interval -> le genereaza ca double * 10^p
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_int_distribution<> generator(copieCapatStanga, capatDreapta);

    return generator(mt);
}

pair<double, int> genereazaValoriDoubleInt(){
    int valoareIndivid = genereazaValoareIntreaga(); // numarul asta il scriem in binar
    int putere = putere10(precizie);

    auto valoareIndividDouble = (double) (valoareIndivid + capatStangaInterval * putere) / putere;

    return {valoareIndividDouble, valoareIndivid};
}

Individ::Individ(double fitnessIndivid, double valoareIndivid, int valoareIntreaga, const vector<int> &cromozom)
        : fitnessIndivid(fitnessIndivid), valoareIndivid(valoareIndivid), valoareIntreaga(valoareIntreaga),
          cromozom(cromozom) {}

vector<Individ> generareIndivizi(){
    vector<Individ> populatie;
    for(int i = 0; i < dimensiuneaPopulatiei; i++){
        auto valori = genereazaValoriDoubleInt();

        double valoareDouble = valori.first;
        int valoareIntreaga = valori.second;

        double fitnessIndivid = fitness(valoareDouble);
        vector<int> cromozom = conversieZecimalInBinar(valoareIntreaga);

        Individ individCurent = Individ(fitnessIndivid, valoareDouble, valoareIntreaga, cromozom);
        populatie.push_back(individCurent);
    }
    return populatie;
}

double sumaFitness(const vector<Individ>& populatie){
    double suma = 0;

    for(const auto& individ : populatie){
        suma += individ.getFitnessIndivid();
    }

    return suma;
}

vector<double> probabilitatiDeSelectie(const vector<Individ>& populatie){
    double sumaFitnessPopulatie = sumaFitness(populatie);
    vector<double> probabilitatiSelectie;

    for(const auto& individ : populatie){
        double probabilitateSelectieIndivid = individ.getFitnessIndivid() / sumaFitnessPopulatie;
        probabilitatiSelectie.push_back(probabilitateSelectieIndivid);
    }

    return probabilitatiSelectie;
}

vector<double> intervaleDeSelectie(vector<double> probabilitatiSelectie){
    vector<double> intervaleSelectie;
    intervaleSelectie.push_back(0);

    for(int i = 0; i < probabilitatiSelectie.size(); i++){
        // elementul curent din intervaleSelectie + elementul de pe pozitia i din probabilitatiSelectie
        intervaleSelectie.push_back(intervaleSelectie[i] + probabilitatiSelectie[i]);
    }

    return intervaleSelectie;
}

Individ determinareElitist(vector<Individ> populatie){
    int pozitieElitist = 0;
    double valoareMaximaFitness = populatie[0].getFitnessIndivid();

    for(int i = 1; i < dimensiuneaPopulatiei; i++){
        if(populatie[i].getFitnessIndivid() > valoareMaximaFitness){
            pozitieElitist = i;
            valoareMaximaFitness = populatie[i].getFitnessIndivid();
        }
    }

    return populatie[pozitieElitist];
}

int cautareBinaraInterval(double u, int capatStanga, int capatDreapta, vector<double> intervaleSelectie){
    if(u >= intervaleSelectie[capatDreapta]){
        return capatDreapta + 1;
    }
    else{
        if(u <= intervaleSelectie[capatStanga]){
            return capatStanga;
        }
        else{
            if(capatStanga < capatDreapta){
                int mijloc = (capatDreapta + capatStanga) / 2;
                if(u >= intervaleSelectie[mijloc] && u < intervaleSelectie[mijloc + 1]){
                    // u apartine intervalului [mijloc, mijloc+1)
                    return mijloc + 1;
                }
                else{
                    if(intervaleSelectie[mijloc] <= u && intervaleSelectie[mijloc + 1] <= u){
                        // u se afla in dreapta intervalului [mijloc, mijloc+1)
                        return cautareBinaraInterval(u, mijloc + 2, capatDreapta, intervaleSelectie);
                    }
                    else{
                        return cautareBinaraInterval(u, capatStanga, mijloc - 1, intervaleSelectie);
                    }
                }
            }
        }
    }
    return 0;
}

// Procesul de selectie
pair<vector<Individ>, vector<pair<int, double>>> selectieCromozomi(vector<Individ> populatie, const vector<double>& intervaleSelectie){
    /*
        - generam un numar aleator u
        - cautam intervalul [qi, q(i+1)), caruia ii apartine u
        - selectam cromozomul q(i+1)
     */

    vector<Individ> indiviziSelectati; // u, indexul, individul
    vector<pair<int, double>> indexIndividSiValoareRandom;

    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> generator(0, 1);

    // selectam individul elitist, deci mai raman n - 1 indivizi de selectat
    for(int i = 0; i < dimensiuneaPopulatiei - 1; i++){
        double u = generator(mt);

        int indexCromozom = cautareBinaraInterval(u, 0, (int)intervaleSelectie.size() - 1, intervaleSelectie);
//        cout << indexCromozom << " ";

        indiviziSelectati.push_back(populatie[indexCromozom - 1]); // vectorul populatie incepe de la 0
        indexIndividSiValoareRandom.emplace_back(indexCromozom - 1, u);
    }

    return {indiviziSelectati, indexIndividSiValoareRandom};
}

vector<pair<bool, double>> selectieIncrucisare(const vector<Individ>& indiviziSelectati){
    /*
        - generam un numar aleator u pentru fiecare individ
        - daca u < probabilitatea de recombinare, cromozomul participa la recombinare
     */

    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> generator(0, 1);

    // indivizi[i] - statutul si valoarea u pentru indiviziSelectati[i]
    vector<pair<bool, double>> indivizi; // returnez un vector de bool (true->se duce la recombinare, false->nu se duce) si double (valoarea u -> retin pentru a o afisa in main)

    for(int i = 0; i < indiviziSelectati.size(); i++){
        double u = generator(mt);

        if (u < probabilitateRecombinare){
            indivizi.emplace_back(true, u);
        }
        else {
            indivizi.emplace_back(false, u);
        }
    }
    return indivizi;
}

int conversieBinarInZecimal(vector<int> numarBinar){
    reverse(numarBinar.begin(), numarBinar.end());
    int numar = 0, putere = 1;
    for(int bit : numarBinar){
        numar += bit * putere;
        putere *= 2;
    }

    return numar;
}

double conversieIntInDouble(int numarIntreg){
    int putere = putere10(precizie);

    auto numarDouble = (double) (numarIntreg + capatStangaInterval * putere) / putere;

    return numarDouble;
}

vector<pair<int, Individ>> recombinare(vector<pair<int, Individ>> indiviziRecombinare, ofstream &out, int indiceEtapa){
    vector<pair<int, Individ>> indiviziDupaRecombinare;

    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_int_distribution<> generator(0, lungimeCromozom - 1);

    for(int i = 0; i < indiviziRecombinare.size(); i = i + 2){
        if(i < indiviziRecombinare.size() && (i + 1) < indiviziRecombinare.size()){
            // in indiviziRecombinare, pe prima pozitie e indicele individului, pe a doua e obiectul de tip individ
            auto cromozomParinte1 = indiviziRecombinare[i].second.getCromozom();
            auto cromozomParinte2 = indiviziRecombinare[i + 1].second.getCromozom();

            vector<int> cromozomCopil1, cromozomCopil2;
            int rupere = generator(mt);

            for(int j = 0; j < rupere; j++){
                cromozomCopil1.push_back(cromozomParinte2[j]);
                cromozomCopil2.push_back(cromozomParinte1[j]);
            }

            for(int j = rupere; j < lungimeCromozom; j++){
                cromozomCopil1.push_back(cromozomParinte1[j]);
                cromozomCopil2.push_back(cromozomParinte2[j]);
            }

            int valoareIntCopil1 = conversieBinarInZecimal(cromozomCopil1), valoareIntCopil2 = conversieBinarInZecimal(cromozomCopil2);
            double valoareDoubleCopil1 = conversieIntInDouble(valoareIntCopil1), valoareDoubleCopil2 = conversieIntInDouble(valoareIntCopil2);

            Individ copil1 = Individ(fitness(valoareDoubleCopil1), valoareDoubleCopil1, valoareIntCopil1, cromozomCopil1);
            Individ copil2 = Individ(fitness(valoareDoubleCopil2), valoareDoubleCopil2, valoareIntCopil2, cromozomCopil2);

            indiviziDupaRecombinare.emplace_back(indiviziRecombinare[i].first, copil1);
            indiviziDupaRecombinare.emplace_back(indiviziRecombinare[i + 1].first, copil2);

            if(indiceEtapa == 0){
                // Doar in prima etapa se fac afisari
                out << "\n - Recombinare dintre cromozomul " << indiviziRecombinare[i].first + 1 << " si cromozomul " << indiviziRecombinare[i + 1].first + 1; // in vector, pozitiile lor incep de la 0
                out << "\n";

                for(int k = 0; k < lungimeCromozom; k++){
                    out << cromozomParinte1[k];
                }

                out << " ";
                for(int k = 0; k < lungimeCromozom; k++){
                    out << cromozomParinte2[k];
                }

                out << " punct de rupere = " << rupere << "\n";
                out << "Rezultat ";

                for(int k = 0; k < lungimeCromozom; k++){
                    out << cromozomCopil1[k];
                }

                out << " ";
                for(int k = 0; k < lungimeCromozom; k++){
                    out << cromozomCopil2[k];
                }
                out << "\n";
            }
        }
        else {
            // in cazul in care lungimea vectorului indiviziRecombinare e impara si ultimul nu mai are cu cine sa se recombine, il adaugam
            indiviziDupaRecombinare.emplace_back(i, indiviziRecombinare[i].second);
        }
    }
    return indiviziDupaRecombinare;
}

vector<pair<int, Individ>> mutatie(vector<Individ> indiviziDupaRecombinare){
    vector<pair<int, Individ>> indiviziMutatie;

    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> generator(0, 1);

    for(int i = 0; i < indiviziDupaRecombinare.size(); i++){
        auto cromozom = indiviziDupaRecombinare[i].getCromozom();
        vector<int> cromozomNou;
        bool mutant = false;

        for(int j = 0; j < cromozom.size(); j++) {
            double u = generator(mt);
            if (u < probabilitateMutatie) {
                if (cromozom[j] == 0) {
                    cromozomNou.push_back(1);
                } else {
                    cromozomNou.push_back(0);
                }
                mutant = true;
            }
            else{
                cromozomNou.push_back(cromozom[j]);
            }
        }

        if (mutant){
            // daca i-am modificat vreo gena, generam un nou individ cu cromozomul nou
            int valoareNouaInt = conversieBinarInZecimal(cromozomNou);
            double valoareNouaDouble = conversieIntInDouble(valoareNouaInt);

            Individ individNou = Individ(fitness(valoareNouaDouble), valoareNouaDouble, valoareNouaInt, cromozomNou);

            indiviziMutatie.emplace_back(i, individNou);
        }
    }

    return indiviziMutatie;
}

/* Functii din clasa Individ */

double Individ::getFitnessIndivid() const {
    return fitnessIndivid;
}

void Individ::setFitnessIndivid(double fitness) {
    Individ::fitnessIndivid = fitness;
}

double Individ::getValoareIndivid() const {
    return valoareIndivid;
}

void Individ::setValoareIndivid(double valoare) {
    Individ::valoareIndivid = valoare;
}

const std::vector<int> &Individ::getCromozom() const {
    return cromozom;
}

void Individ::setCromozom(const std::vector<int> &cromozomIndivid) {
    Individ::cromozom = cromozomIndivid;
}

int Individ::getValoareIntreaga() const {
    return valoareIntreaga;
}

void Individ::setValoareIntreaga(int valoareIntreagaIndivid) {
    Individ::valoareIntreaga = valoareIntreagaIndivid;
}

void afisareIndivid(const Individ& individ, ofstream &out, int indice){
    out << indice + 1 << ": ";
    auto cromozomCurent = individ.getCromozom();

    for(int k = 0; k < lungimeCromozom; k++){
        out << cromozomCurent[k];
    }

    out << " x = " << individ.getValoareIndivid();
    out << " f = " << individ.getFitnessIndivid() << "\n";
}

/*
 In fisier, datele de intrare sunt de forma:
  - dimensiunea populatiei
  - domeniul de definitie al functiei
  - coeficientii polinomului de gradul 2
  - precizia
  - probabilitatea de recombinare
  - probabilitatea de mutatie
  - numarul de etape al algoritmului
 */

int main() {
    ifstream in("problema.in");
    ofstream out("Evolutie.out");

    in >> dimensiuneaPopulatiei >> capatStangaInterval >> capatDreaptaInterval >> a >> b >> c >> precizie >> probabilitateRecombinare >> probabilitateMutatie >> numarEtape;
    lungimeCromozom = calculLungimeCromozom();

    vector<Individ> populatie = generareIndivizi();
//    vector<Individ> elitisti;
    vector<double> valoriMediiPerformanta;

    for(int i = 0; i < numarEtape; i++){
        vector<double> probabilitatiSelectie = probabilitatiDeSelectie(populatie);
        vector<double> intervaleSelectie = intervaleDeSelectie(probabilitatiSelectie);

        auto selectie = selectieCromozomi(populatie, intervaleSelectie);
        vector<Individ> indiviziSelectati = selectie.first;
        vector<pair<int, double>> valori = selectie.second;

        auto valoriIncrucisare = selectieIncrucisare(indiviziSelectati);
        vector<pair<int, Individ>> indiviziRecombinare; // index si individ

        for (int j = 0; j < valoriIncrucisare.size(); j++) {
            if (valoriIncrucisare[j].first) {
                indiviziRecombinare.emplace_back(j, indiviziSelectati[j]);
            }
        }

        if(i == 0) {
            out << " * Populatia initiala\n";
            for (int j = 0; j < dimensiuneaPopulatiei; j++) {
                afisareIndivid(populatie[j], out, j);
            }

            out << "\n * Probabilitati selectie\n";
            for (int j = 0; j < dimensiuneaPopulatiei; j++) {
                out << "cromozom " << j + 1;
                out << " probabilitate " << probabilitatiSelectie[j] << "\n";
            }

            out << "\n * Intervale probabilitati selectie\n";
            for (auto interval: intervaleSelectie) {
                out << interval << " ";
            }
            out << "\n";

            for (int j = 0; j < indiviziSelectati.size(); j++) {
                out << "\nu = " << valori[j].second << " selectam cromozomul " << valori[j].first + 1;
            }

            out << "\n\n * Dupa selectie:\n";

            for (int j = 0; j < dimensiuneaPopulatiei - 1; j++) {
                afisareIndivid(indiviziSelectati[j], out, j);
            }

            out << "\n * Probabilitatea de incrucisare = " << probabilitateRecombinare << "\n";

            for (int j = 0; j < valoriIncrucisare.size(); j++) {
                out << j + 1 << ": ";
                auto cromozomCurent = indiviziSelectati[j].getCromozom();

                for (int k = 0; k < lungimeCromozom; k++) {
                    out << cromozomCurent[k];
                }

                out << " u = " << valoriIncrucisare[j].second;
                if (valoriIncrucisare[j].first) {
                    out << " < " << probabilitateRecombinare << " participa la recombinare";
                }
                out << "\n";
            }
        }

        // functia recombinare afiseaza si indivizii recombinati
        auto indiviziRecombinati = recombinare(indiviziRecombinare, out, i);
        auto indiviziDupaRecombinare = indiviziSelectati;

        for(auto &individ : indiviziRecombinati){
            int indiceRecombinat = individ.first;
            indiviziDupaRecombinare[indiceRecombinat] = individ.second;
        }

        // indivizii dupa recombinare sunt stocati in vectorul indiviziSelectati!
        auto indiviziDupaMutatie = mutatie(indiviziDupaRecombinare);
        auto indiviziMutanti = indiviziDupaRecombinare;

        for(auto &individ : indiviziDupaMutatie){
            int indiceMutant = individ.first;
            indiviziMutanti[indiceMutant] = individ.second;
        }

        if(i == 0){
            out << "\n * Dupa etapa de recombinare: \n";
            for(int j = 0; j < indiviziDupaRecombinare.size(); j++){
                afisareIndivid(indiviziDupaRecombinare[j], out, j);
            }

            out << "\n * Probabilitatea de mutatie pentru fiecare gena = " << probabilitateMutatie << "\n";
            for(auto &individ : indiviziDupaMutatie){
                int indiceMutant = individ.first;
                out << indiceMutant + 1 << "\n";
            }

            out << " * Dupa mutatie: \n";
            for(int j = 0; j < indiviziMutanti.size(); j++){
                afisareIndivid(indiviziMutanti[j], out, j);
            }

            out << "\n\n * Evolutia maximului\n";
        }

        Individ individElitistCurent = determinareElitist(populatie);
//        elitisti.push_back(individElitistCurent);

        vector<Individ> generatieNoua = indiviziMutanti;
        generatieNoua.push_back(individElitistCurent);

        populatie = generatieNoua;

        out << setprecision(25) << individElitistCurent.getFitnessIndivid() << "\n";

        double valoareMediePerformanta = sumaFitness(populatie) / dimensiuneaPopulatiei;
        valoriMediiPerformanta.push_back(valoareMediePerformanta);
    }

    out << "\n * Valori medii performanta (pentru fiecare etapa)\n";
    for(int i = 0; i < numarEtape; i++){
        out << "Etapa " << i + 1 << ": " << valoriMediiPerformanta[i] << "\n";
    }

//    cout << generator(mt) << " " << generator(mt) << " " << generator(mt);

    in.close();
    out.close();
    return 0;
}
