#include <filesystem>         // Pour parcourir les fichiers dans un dossier
#include <chrono>             // Pour mesurer le temps d'exécution
#include <fstream>            // Pour lire/écrire des fichiers
#include <iostream>           // Pour afficher des messages à la console
#include <vector>             // Pour les vecteurs dynamiques
#include <random>             // Pour générer des valeurs aléatoires
#include <ctime>              // Pour utiliser l'heure comme graine (seed)
#include <nlohmann/json.hpp>  // Pour parser les fichiers JSON
#include <set>                // Pour des ensembles de valeurs uniques
#include <iomanip>            // Pour la mise en forme des sorties
#include <numeric>            // Pour accumulate, etc.
#include <cmath>              // Pour les fonctions mathématiques


using namespace std;
using json = nlohmann::json;

using Chromosome = vector<int>;
using Population = vector<Chromosome>;

int N_CAPTEURS = 0;
int N_EMPLACEMENTS = 0;
int TAILLE_POPULATION = 0;
int TAILLE_CHROMOSOME = N_CAPTEURS * N_EMPLACEMENTS;
int MAX_GENERATIONS = 3;

vector<pair<double, double>> emplacements;         // Coordonnées des emplacements possibles
vector<pair<double, double>> points_interet;       // Coordonnées des POI à couvrir
vector<vector<double>> matrice_distance;           // Distances entre emplacements et POI
vector<double> rayons_capteurs;                    // Rayon de chaque capteur


// Lecture des données JSON
void lireDonneesDepuisJSON(const string& nomFichier) {
    ifstream fichier(nomFichier);
    if (!fichier) {
        cerr << "Erreur ouverture du fichier JSON.\n";// Affiche une erreur si le fichier est introuvable
        exit(1);
    }
    json data;
    fichier >> data;// Parse le contenu du JSON
 // Lecture des paramètres
    N_CAPTEURS = data["capteurs"].size();
    N_EMPLACEMENTS = data["emplacements"].size();
    TAILLE_POPULATION = data["genetic_params"]["population_size"];
    MAX_GENERATIONS = data["genetic_params"]["generations"];
    TAILLE_CHROMOSOME = N_CAPTEURS;

    emplacements.clear();
    points_interet.clear();
    rayons_capteurs.clear();
    // Chargement des positions et rayon
    for (const auto& e : data["emplacements"])
        emplacements.emplace_back(e["x"], e["y"]);

    for (const auto& p : data["points_interet"])
        points_interet.emplace_back(p["x"], p["y"]);

    for (const auto& capteur : data["capteurs"])
        rayons_capteurs.push_back(capteur["rayon"]);
}

// Calcul de la distance euclidienne
double calculerDistance(const pair<double, double>& a, const pair<double, double>& b) {
    return sqrt(pow(a.first - b.first, 2) + pow(a.second - b.second, 2));
}

// Calcul de la matrice de distances
void calculerMatriceDistance() {
    matrice_distance.resize(N_EMPLACEMENTS, vector<double>(points_interet.size()));
    for (int i = 0; i < N_EMPLACEMENTS; ++i)
        for (int j = 0; j < points_interet.size(); ++j)
            matrice_distance[i][j] = calculerDistance(emplacements[i], points_interet[j]);
}
//Cette fonction vérifie si tous les points d’intérêt ont bien été couverts par les capteurs dans une solution donnée (chromosome).
//C’est une condition essentielle pour qu’une solution soit considérée comme "valide".
// Vérifie si la solution couvre tous les POI
bool estCouvertureComplete(const Chromosome& ind) {
    vector<bool> couverts(points_interet.size(), false);//On crée une liste de "vrai/faux" pour chaque point d'intérêt (POI), tous initialisés à false, car au départ, aucun point n'est encore couvert.
    for (int i = 0; i < N_CAPTEURS; ++i) {
        if (ind[i] == 0) continue;//Si ce capteur n’est pas utilisé dans cette solution (valeur 0), on le saute.
        int e = ind[i] - i * N_EMPLACEMENTS - 1;//calcule lemplacement de capteur utiliser et il retounre l'indice  de l’emplacement
        if (e < 0 || e >= matrice_distance.size()) continue;
        double rayon = rayons_capteurs[i];//On récupère le rayon de couverture du capteur
        for (int j = 0; j < points_interet.size(); ++j)
            if (!couverts[j] && matrice_distance[e][j] <= rayon)//Si ce point n'est pas encore couvert et qu’il est dans le rayon du capteur → alors on le marque comme couvert.
                couverts[j] = true;
    }
    return all_of(couverts.begin(), couverts.end(), [](bool b) { return b; });//vérifie si tous les points sont couverts.  si oui rentur yes sinon f
}
//Cette fonction récompense les solutions qui : Couvre tous les points d’intérêt Utilisent le moins de capteurs possible 
// Fonction d'évaluation
double fonctionFitness(const Chromosome& ind) {//Elle retourne une valeur double qui représente la performance de ce chromosome
    if (!estCouvertureComplete(ind)) return 0.0;//Si la solution ne couvre pas tous les points d’intérêt, on retourne immédiatement 0.0.
//Ça évite de favoriser des solutions incomplètes — c’est un filtre pour ne garder que les solutions valides.

    int capteurs_utiles = 0;// On initialise un compteur pour savoir combien de capteurs ont été réellement utilisés dans cette solution.
   // On compte les capteurs activés (valeurs ≠ 0) pour savoir combien sont utilisés dans la solution, afin de pouvoir ensuite appliquer une pénalité si leur nombre est élevé
    for (int i = 0; i < N_CAPTEURS; ++i)
        if (ind[i] != 0) capteurs_utiles++;

    double penalite = 0.8 * capteurs_utiles;//On diminue le score (fitness) en fonction du nombre de capteurs utilisés — plus on utilise de capteurs, plus on est pénalisé.
   //Le nombre total de points d’intérêt (objectif à maximiser) moins la pénalité (objectif à minimiser)  || Cela équilibre la couverture avec l’économie de ressources (capteurs).


    return static_cast<double>(points_interet.size()) - penalite;
}

// Initialisation de la population
Population initialiserPopulation() {
    Population population;
    random_device rd;
    mt19937 gen(rd());

    for (int i = 0; i < TAILLE_POPULATION; ++i) {
        Chromosome individu;
        for (int j = 0; j < N_CAPTEURS; ++j) {
            uniform_real_distribution<> chance(0.0, 1.0);
            if (chance(gen) < 0.3)
                individu.push_back(0);
            else {
                int base = j * N_EMPLACEMENTS;
                uniform_int_distribution<> distrib(base + 1, base + N_EMPLACEMENTS);
                individu.push_back(distrib(gen));
            }
        }
        population.push_back(individu);
    }
    return population;
}

// Sélection par tournoi
Chromosome tournoiSelection(const Population& population, int k = 5) {
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> distrib(0, population.size() - 1);
    Chromosome meilleur;
    double meilleure_fitness = -1e9;
    for (int i = 0; i < k; ++i) {
        int idx = distrib(gen);
        double f = fonctionFitness(population[idx]);
        if (f > meilleure_fitness) {
            meilleure_fitness = f;
            meilleur = population[idx];
        }
    }
    return meilleur;
}

// Croisement à un point
pair<Chromosome, Chromosome> croisement(const Chromosome& p1, const Chromosome& p2) {
    int taille = p1.size();
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dist(1, taille - 2);
    int point = dist(gen);
    Chromosome e1 = p1, e2 = p2;
    for (int i = point; i < taille; ++i)
        swap(e1[i], e2[i]);
    return {e1, e2};
}

// Mutation
void mutation(Chromosome& ind, double taux = 0.1) {
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> p(0.0, 1.0);
    for (int i = 0; i < ind.size(); ++i) {
        if (p(gen) < taux) {
            int base = i * N_EMPLACEMENTS;
            uniform_int_distribution<> d(base + 1, base + N_EMPLACEMENTS);
            ind[i] = d(gen);
        }
    }
}

// Écriture des résultats CSV
void enregistrerResultat(const string& nom, int nb_poi, int nb_caps, bool faisable,
                         int utilises, int nb_couverts, long long temps, int gen, const string& csv) {
    ofstream f(csv, ios::app);
    f << nom << "," << nb_poi << "," << nb_caps << "," << (faisable ? "Oui" : "Non") << ","
      << utilises << "," << nb_couverts << "," << temps << "," << gen << "\n"; // ajout nb_couverts
    f.close();
}


// Initialisation fichier CSV
void initialiserFichierResultats(const string& nomFichier) {
    ofstream f(nomFichier);
    f << "Fichier,POI,Capteurs,Faisable,Capteurs_utilisés,POI_Couverts,Temps_ms,Generations\n"; // ajouté POI_Couverts
    f.close();
}


// MAIN
int main() {
//Définit le chemin du dossier contenant les fichiers d'entrée JSON et le nom du fichier CSV de sortie
    const string dossier = "experiences/";
    const string fichier_csv = "resultats10.csv";
//Vérifie si le fichier CSV existe. Sinon, appelle initialiserFichierResultats() pour créer le fichier et écrire l'en-tête (supposé contenir des colonnes comme nom_fichier, faisable
    if (!filesystem::exists(fichier_csv))
        initialiserFichierResultats(fichier_csv);
//Itère sur tous les fichiers du dossier experiences/ pour traiter chaque expérience
    for (const auto& entry : filesystem::directory_iterator(dossier)) {
 //Récupère le chemin complet du fichier et enregistre l'heure de début pour mesurer le temps d'exécution.   
        string nom_fichier = entry.path().string();
        auto debut = chrono::high_resolution_clock::now();
//Nettoie les variables globales (supposées) avant de traiter un nouveau fichier pour éviter des données résiduelles.
        emplacements.clear();
        points_interet.clear();
        rayons_capteurs.clear();
        matrice_distance.clear();
//Lit les emplacements, points d'intérêt et rayons des capteurs depuis le fichier JSON.
        lireDonneesDepuisJSON(nom_fichier);
        calculerMatriceDistance();//Calcule la matrice de distance entre chaque emplacement et chaque point d'intérêt (pour optimiser les calculs ultérieurs).
        Population population = initialiserPopulation();
// Initialise les variables pour stocker la meilleure solution trouvée, sa fitness, si elle est faisable (couverture complète), et la génération où elle a été trouvée.
        Chromosome meilleur;
        double best_fitness = -1e9;
        bool faisable = false;
        int generation_finale = -1;
//Boucle sur un nombre fixe de générations (MAX_GENERATIONS).
        for (int gen = 1; gen <= MAX_GENERATIONS; ++gen) {
//    Évalue chaque individu avec fonctionFitness(). //Met à jour la meilleure solution si une fitness supérieure est trouvée. //Vérifie si la solution est faisable via estCouvertureComplete() (couverture totale des points d'intérêt).          
            for (const Chromosome& ind : population) {
                double fit = fonctionFitness(ind);
                if (fit > best_fitness) {
                    best_fitness = fit;
                    meilleur = ind;
                    generation_finale = gen;
                    faisable = estCouvertureComplete(ind);
                }
            }

            if (faisable) break;//rrête l'algorithme si une solution faisable est trouvée.
//Sélectionne des individus par tournoi pour la reproduction.
            Population nouvelle;
            for (int i = 0; i < TAILLE_POPULATION; ++i)
                nouvelle.push_back(tournoiSelection(population));
//Applique un opérateur de croisement pour générer de nouveaux individus.
            Population croisee;
            for (int i = 0; i < TAILLE_POPULATION; i += 2) {
                auto [e1, e2] = croisement(nouvelle[i], nouvelle[(i + 1) % TAILLE_POPULATION]);
                croisee.push_back(e1);
                croisee.push_back(e2);
            }
 //Applique une mutation aléatoire et remplace l'ancienne population.	           
            for (Chromosome& ind : croisee)
                mutation(ind);
            population = croisee;
        }
//Calcule le temps d'exécution en millisecondes.
        auto fin = chrono::high_resolution_clock::now();
        auto temps_ms = chrono::duration_cast<chrono::milliseconds>(fin - debut).count();
//Compte le nombre de capteurs actifs dans la meilleure solution.
        int capteurs_utiles = 0;
        for (int i = 0; i < N_CAPTEURS; ++i)
            if (meilleur[i] != 0) capteurs_utiles++;
//Si la solution est faisable, calcule le nombre de points d'intérêt couverts. Utilise la matrice_distance pour vérifier la couverture par chaque capteur actif.            
            int nb_couverts = 0;
if (faisable) {
    vector<bool> couverts(points_interet.size(), false);
    for (int i = 0; i < N_CAPTEURS; ++i) {
        if (meilleur[i] == 0) continue;
        int e = meilleur[i] - i * N_EMPLACEMENTS - 1;
        double rayon = rayons_capteurs[i];
        for (int j = 0; j < points_interet.size(); ++j) {
            if (!couverts[j] && matrice_distance[e][j] <= rayon)
                couverts[j] = true;
        }
    }
    nb_couverts = count(couverts.begin(), couverts.end(), true);
}
//Enregistre les résultats dans le CSV (nom du fichier, nombre de points, capteurs, faisabilité
        enregistrerResultat(
    entry.path().filename().string(),
    points_interet.size(),
    N_CAPTEURS,
    faisable,
    capteurs_utiles,
    nb_couverts,
    temps_ms,
    generation_finale == -1 ? MAX_GENERATIONS : generation_finale,
    fichier_csv
);

    }
// Affiche un message de succès et termine le programm
    cout << "\n Toutes les expériences ont été traitées.\n";
    return 0;
}
