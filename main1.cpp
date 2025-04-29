#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <set>
#include <random>
#include "json.hpp"

using namespace std;
using json = nlohmann::json;

// Structures
struct Position {
    double x;
    double y;
};

struct PointInteret {
    int id;
    Position pos;
};

struct Capteur {
    int id;
    double rayon;
};

// Variables globales
int TAILLE_POPULATION;
int N_GENERATIONS;
double MUTATION_RATE;
double CROSSOVER_RATE;
int TOURNAMENT_SIZE;
int N_CAPTEURS = 0;
int N_EMPLACEMENTS = 0;
int chromosome_size = 0;

vector<Capteur> capteurs;
vector<PointInteret> points_interet;
vector<Position> emplacements;
vector<vector<double>> matrice_distance;

// Types
typedef vector<int> Chromosome;
typedef vector<Chromosome> Population;
void reparer(Chromosome& chromo);
// Logger
void log(string texte) {
    ofstream fichier("log.txt", ios::app);
    if (fichier.is_open()) {
        fichier << texte << endl;
        fichier.close();
    }
}

// Fonctions
double calculerDistance(Position p1, Position p2) {
    return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2));
}

void calculerMatriceDistance() {
    matrice_distance.resize(emplacements.size(), vector<double>(points_interet.size()));
    for (size_t i = 0; i < emplacements.size(); ++i) {
        for (size_t j = 0; j < points_interet.size(); ++j) {
            matrice_distance[i][j] = calculerDistance(emplacements[i], points_interet[j].pos);
        }
    }

    log("Matrice de distance calculée :");
    for (size_t i = 0; i < emplacements.size(); ++i) {
        string ligne = "Emplacement " + to_string(i+1) + ": ";
        for (size_t j = 0; j < points_interet.size(); ++j) {
            ligne += to_string(matrice_distance[i][j]) + " ";
        }
        log(ligne);
    }
}

void calculerEtAfficherMatriceCouvertureParCapteur() {
    for (size_t capteur_idx = 0; capteur_idx < capteurs.size(); ++capteur_idx) {
        double rayon = capteurs[capteur_idx].rayon;

        vector<vector<int>> matrice_couverture(emplacements.size(), vector<int>(points_interet.size(), 0));

        for (size_t i = 0; i < emplacements.size(); ++i) {
            for (size_t j = 0; j < points_interet.size(); ++j) {
                if (matrice_distance[i][j] <= rayon) {
                    matrice_couverture[i][j] = 1;
                }
            }
        }

        // Afficher cette matrice
        log("=== Matrice de couverture du Capteur " + to_string(capteurs[capteur_idx].id) + " (rayon = " + to_string(rayon) + ") ===");
        for (size_t i = 0; i < emplacements.size(); ++i) {
            string ligne = "Emp " + to_string(i+1) + ": ";
            for (size_t j = 0; j < points_interet.size(); ++j) {
                ligne += to_string(matrice_couverture[i][j]) + " ";
            }
            log(ligne);
        }
        log("===============================================================");
    }
}

vector<vector<int>> initialiser_population() {
    vector<vector<int>> population;
    random_device rd;
    mt19937 gen(rd());

    int nb_emplacements_par_capteur = N_EMPLACEMENTS;
    for (int i = 0; i < TAILLE_POPULATION; i++) {
        vector<int> chromosome(chromosome_size, 0);

        for (int capteur = 0; capteur < N_CAPTEURS; capteur++) {
            int debut = capteur * nb_emplacements_par_capteur;
            uniform_int_distribution<> dis(debut, debut + nb_emplacements_par_capteur - 1);
            int position = dis(gen);
            chromosome[position] = 1; // activer un seul bit par capteur
        }

        population.push_back(chromosome);
    }

    return population;
}

double fonctionFitness(const Chromosome& individu);
void afficher_population(const vector<vector<int>>& population) {
    ofstream fichier("log.txt", ios::app); // ouvrir en mode ajout

    log("Population initiale :");
    for (int i = 0; i < population.size(); i++) {
        string ligne = "Individu " + to_string(i + 1) + " : ";
        bool first = true;
        for (int j = 0; j < population[i].size(); j++) {
            if (population[i][j] == 1) {
                if (!first) ligne += ", ";
                ligne += to_string(j + 1); // commencer à 1
                first = false;
            }
        }

        double fitness = fonctionFitness(population[i]);
        ligne += " | Fitness = " + to_string(fitness);

        cout << ligne << endl;         // afficher sur la console
        if (fichier.is_open()) {
            fichier << ligne << endl;  // écrire aussi dans log.txt
        }
    }

    log("--------------------------------------------------");
    fichier.close(); // fermer le fichier à la fin
}





bool estCouvertureComplete(const Chromosome& individu) {
    set<int> couverts;

    int nb_emplacements_par_capteur = N_EMPLACEMENTS;
    for (int capteur = 0; capteur < N_CAPTEURS; capteur++) {
        int debut = capteur * nb_emplacements_par_capteur;
        int id_emplacement = -1;

        for (int pos = debut; pos < debut + nb_emplacements_par_capteur; ++pos) {
            if (individu[pos] == 1) {
                id_emplacement = pos - debut; // Position locale (0..N-1)
                break;
            }
        }

        if (id_emplacement == -1) continue; // Aucun capteur placé

        for (size_t j = 0; j < points_interet.size(); ++j) {
            if (matrice_distance[id_emplacement][j] <= capteurs[capteur].rayon) {
                couverts.insert(j);
            }
        }
    }
    return couverts.size() == points_interet.size();
}/*
bool estConnecte(const Chromosome& individu) {
    vector<pair<Position, double>> capteurs_actifs; // (position, rayon)

    int nb_emplacements_par_capteur = N_EMPLACEMENTS;

    for (int capteur = 0; capteur < N_CAPTEURS; ++capteur) {
        double rayon = capteurs[capteur].rayon;
        if (rayon == 0) continue; // ignorer capteur inutile

        int debut = capteur * nb_emplacements_par_capteur;
        for (int pos = debut; pos < debut + nb_emplacements_par_capteur; ++pos) {
            if (individu[pos] == 1) {
                capteurs_actifs.push_back({emplacements[pos - debut], rayon});
                break;
            }
        }
    }

    if (capteurs_actifs.empty()) return false; // Aucun capteur actif

    // Construire le graphe
    int n = capteurs_actifs.size();
    vector<vector<int>> graph(n);

    for (int i = 0; i < n; ++i) {
        Position p1 = capteurs_actifs[i].first;
        double r1 = capteurs_actifs[i].second;
        for (int j = i + 1; j < n; ++j) {
            Position p2 = capteurs_actifs[j].first;
            double r2 = capteurs_actifs[j].second;
            double dist = calculerDistance(p1, p2);
            if (dist <= r1 + r2) {
                graph[i].push_back(j);
                graph[j].push_back(i);
            }
        }
    }

    // Parcours DFS pour vérifier la connectivité
    vector<bool> visité(n, false);
    vector<int> stack;
    stack.push_back(0);
    visité[0] = true;

    while (!stack.empty()) {
        int u = stack.back();
        stack.pop_back();
        for (int v : graph[u]) {
            if (!visité[v]) {
                visité[v] = true;
                stack.push_back(v);
            }
        }
    }

    // Tous les capteurs doivent être visités
    for (bool v : visité) {
        if (!v) return false;
    }

    return true;
}
*/
double fonctionFitness(const Chromosome& individu) {
    if (!estCouvertureComplete(individu)) return 0.0;

    // Connexité désactivée ici

    vector<bool> couverts(points_interet.size(), false);
    int capteurs_utiles = 0;

    int nb_emplacements_par_capteur = N_EMPLACEMENTS;

    for (int capteur = 0; capteur < N_CAPTEURS; ++capteur) {
        double rayon = capteurs[capteur].rayon;
        if (rayon == 0) continue;

        int debut = capteur * nb_emplacements_par_capteur;
        int id_emplacement = -1;

        for (int pos = debut; pos < debut + nb_emplacements_par_capteur; ++pos) {
            if (individu[pos] == 1) {
                id_emplacement = pos - debut;
                break;
            }
        }

        if (id_emplacement == -1) continue;

        bool a_couvert = false;
        for (size_t j = 0; j < points_interet.size(); ++j) {
            if (!couverts[j] && matrice_distance[id_emplacement][j] <= rayon) {
                couverts[j] = true;
                a_couvert = true;
            }
        }

        if (a_couvert) capteurs_utiles++;
    }

    int nb_couverts = count(couverts.begin(), couverts.end(), true);
    double penalite = 0.8 * capteurs_utiles;

    return static_cast<double>(nb_couverts) - penalite;
}


void lireJSON(string nom_fichier) {
    ifstream input(nom_fichier);
    json j;
    input >> j;

    for (auto& c : j["capteurs"]) {
        if (c["rayon"] > 0) {  // <<< ajoute cette condition
            capteurs.push_back({c["id"], c["rayon"]});
        } else {
            log(" Capteur ID " + to_string(c["id"]) + " ignoré (rayon = 0).");
        }
    }

    for (auto& p : j["points_interet"]) {
        points_interet.push_back({p["id"], {p["x"], p["y"]}});
    }

    for (auto& e : j["emplacements"]) {
        emplacements.push_back({e["x"], e["y"]});
    }

    TAILLE_POPULATION = j["genetic_params"]["population_size"];
    N_GENERATIONS = j["genetic_params"]["generations"];
    MUTATION_RATE = j["genetic_params"]["mutation_rate"];
    CROSSOVER_RATE = j["genetic_params"]["crossover_rate"];
    TOURNAMENT_SIZE = j["genetic_params"]["tournament_size"];

    N_CAPTEURS = capteurs.size();
    N_EMPLACEMENTS = emplacements.size();
    chromosome_size = N_CAPTEURS * N_EMPLACEMENTS;
}


//
Chromosome tournoi_selection(const Population& population) {
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(0, population.size() - 1);

    Chromosome meilleur;
    int meilleur_index = -1;
    double meilleure_fitness = -1.0;
log("--------------------------------------------------");
    log("Début de la sélection par tournoi:");
    for (int i = 0; i < TOURNAMENT_SIZE; i++) {
        int index = dis(gen);
        const Chromosome& candidat = population[index];
        double fitness = fonctionFitness(candidat);

        string ligne = "Candidat " + to_string(i+1) + " (Individu " + to_string(index+1) + ") : fitness = " + to_string(fitness);
        log(ligne);

        if (fitness > meilleure_fitness) {
            meilleur = candidat;
            meilleure_fitness = fitness;
            meilleur_index = index;
        }
    }
    log("Fin de la sélection par tournoi.");
    
    if (meilleur_index != -1) {
        log("Individu sélectionné : Individu " + to_string(meilleur_index + 1) + " avec fitness = " + to_string(meilleure_fitness));
    } else {
        log("Erreur : Aucun individu sélectionné !");
    }

    log("--------------------------------------------------");

    return meilleur;
}


tuple<Chromosome, Chromosome, int> croisement(const Chromosome& parent1, const Chromosome& parent2) {
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(0, chromosome_size - 1);

    int point = dis(gen);

    Chromosome enfant1 = parent1;
    Chromosome enfant2 = parent2;

    for (int i = point; i < chromosome_size; i++) {
        swap(enfant1[i], enfant2[i]);
        
    }
reparer(enfant1); // <-- ajout !
reparer(enfant2); // <-- ajout !

    // Construire le texte à logger et afficher
    stringstream ss;
    ss << "\n=== Croisement effectué ===\n";
    ss << "Point de croisement : " << point << "\n";

    ss << "Parent 1 : ";
    for (int gene : parent1) ss << gene << " ";
    ss << "\n";

    ss << "Parent 2 : ";
    for (int gene : parent2) ss << gene << " ";
    ss << "\n";

    ss << "Enfant 1 : ";
    for (int gene : enfant1) ss << gene << " ";
    ss << "\n";

    ss << "Enfant 2 : ";
    for (int gene : enfant2) ss << gene << " ";
    ss << "\n";

    ss << "============================\n";

    // Afficher dans la console
    cout << ss.str() << endl;

    // Sauvegarder dans le fichier
    log(ss.str());
    // Après mutation ou croisement
// Vérification enfant1
for (int capteur = 0; capteur < N_CAPTEURS; capteur++) {
    int nb1 = 0;
    int debut = capteur * N_EMPLACEMENTS;
    for (int j = debut; j < debut + N_EMPLACEMENTS; j++) {
        if (enfant1[j] == 1) nb1++;
    }
    if (nb1 != 1) {
        log("ERREUR dans enfant1 : Le capteur " + to_string(capteur) + " a " + to_string(nb1) + " emplacements actifs !");
    }
}

// Vérification enfant2
for (int capteur = 0; capteur < N_CAPTEURS; capteur++) {
    int nb1 = 0;
    int debut = capteur * N_EMPLACEMENTS;
    for (int j = debut; j < debut + N_EMPLACEMENTS; j++) {
        if (enfant2[j] == 1) nb1++;
    }
    if (nb1 != 1) {
        log(" ERREUR dans enfant2 : Le capteur " + to_string(capteur) + " a " + to_string(nb1) + " emplacements actifs !");
    }
}


    return {enfant1, enfant2, point};
    
}





//
void mutation(Chromosome& individu, int numero_individu) {
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0.0, 1.0);

    Chromosome avant_mutation = individu; // Copier l'individu avant mutation
    vector<int> points_mutation;          // Liste des indices mutés

    for (int i = 0; i < chromosome_size; i++) {
        if (dis(gen) < MUTATION_RATE) {
            individu[i] = 1 - individu[i]; // Inverser le bit
            points_mutation.push_back(i);  // Stocker le point de mutation
        }
    }
    reparer(individu);
    
// Après mutation ou croisement
for (int capteur = 0; capteur < N_CAPTEURS; capteur++) {
    int nb1 = 0;
    int debut = capteur * N_EMPLACEMENTS;
    for (int j = debut; j < debut + N_EMPLACEMENTS; j++) {
        if (individu[j] == 1) nb1++;
    }
    if (nb1 != 1) {
        log(" ERREUR : Le capteur " + to_string(capteur) + " a " + to_string(nb1) + " emplacements actifs !");
    }
}

    stringstream ss; // ➔ UNE SEULE déclaration ici
    ss << "\n=== Mutation de l'individu " << numero_individu << " ===\n";
    ss << "Individu avant mutation : ";
    for (int gene : avant_mutation) ss << gene << " ";
    ss << "\n";

    if (!points_mutation.empty()) {
        ss << "Points de mutation (indices modifiés) : ";
        for (int idx : points_mutation) ss << idx << " ";
        ss << "\n";
    } else {
        ss << "Aucune mutation effectuée.\n";
    }

    ss << "Individu après mutation : ";
    for (int gene : individu) ss << gene << " ";
    ss << "\n=================\n";

    cout << ss.str();   // Afficher sur console
    log(ss.str());      // Sauvegarder dans log.txt
     ss << "\n=================\n";
}

//
Population nouvelle_generation(const Population& population) {
    Population nouvelle_pop;

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0.0, 1.0);

    while (nouvelle_pop.size() < TAILLE_POPULATION) {
        Chromosome parent1 = tournoi_selection(population);
        Chromosome parent2 = tournoi_selection(population);

        Chromosome enfant1 = parent1;
        Chromosome enfant2 = parent2;

        if (dis(gen) < CROSSOVER_RATE) {
            tie(enfant1, enfant2, ignore) = croisement(parent1, parent2);
        }

        mutation(enfant1, nouvelle_pop.size() + 1);
reparer(enfant1);
nouvelle_pop.push_back(enfant1);

if (nouvelle_pop.size() < TAILLE_POPULATION) {
    mutation(enfant2, nouvelle_pop.size() + 1);
    reparer(enfant2);
    nouvelle_pop.push_back(enfant2);
}

    }

    log("Nouvelle génération créée.");
    log("--------------------------------------------------");

    return nouvelle_pop;
}


//
void afficher_chromosome(const Chromosome& chromo) {
    for (bool gene : chromo) {
        std::cout << gene << " ";
    }
    std::cout << std::endl;
}
//

//////////////////////////////////////////////////////////////////////////////

void eliminer_capteurs_inferieurs(Chromosome& chromo) {
    for (int e = 0; e < N_EMPLACEMENTS; ++e) {
        int meilleur_capteur = -1;
        double meilleur_rayon = -1.0;

        // Chercher le capteur actif avec le plus grand rayon pour cet emplacement
        for (int capteur = 0; capteur < N_CAPTEURS; ++capteur) {
            int idx = capteur * N_EMPLACEMENTS + e;
            if (chromo[idx] == 1 && capteurs[capteur].rayon > meilleur_rayon) {
                meilleur_capteur = capteur;
                meilleur_rayon = capteurs[capteur].rayon;
            }
        }

        // Supprimer les autres capteurs de cet emplacement
        for (int capteur = 0; capteur < N_CAPTEURS; ++capteur) {
            int idx = capteur * N_EMPLACEMENTS + e;
            if (capteur != meilleur_capteur) {
                chromo[idx] = 0;
            }
        }
    }
    eliminer_capteurs_inferieurs(chromo);
}


void reparer(Chromosome& chromo) {
    int nb_emplacements_par_capteur = N_EMPLACEMENTS;

    for (int capteur = 0; capteur < N_CAPTEURS; ++capteur) {
        int debut = capteur * nb_emplacements_par_capteur;
        int fin = debut + nb_emplacements_par_capteur;

        int compte_1 = 0;
        vector<int> positions_1;

        for (int i = debut; i < fin; ++i) {
            if (chromo[i] == 1) {
                compte_1++;
                positions_1.push_back(i);
            }
        }

        if (compte_1 == 0) {
            // Aucun emplacement actif ➔ choisir un emplacement aléatoire
            random_device rd;
            mt19937 gen(rd());
            uniform_int_distribution<> dis(debut, fin - 1);
            int idx = dis(gen);
            chromo[idx] = 1;
        } else if (compte_1 > 1) {
            // Plus d'un 1 ➔ en garder 1 au hasard, mettre 0 aux autres
            random_device rd;
            mt19937 gen(rd());
            shuffle(positions_1.begin(), positions_1.end(), gen);
            for (size_t k = 1; k < positions_1.size(); ++k) {
                chromo[positions_1[k]] = 0;
            }
        }
    }
}
//bool estSolutionValide(const Chromosome& individu) {
   // return estCouvertureComplete(individu) && estConnecte(individu);
//}
//
// Écriture des résultats CSV
void enregistrerResultat(const string& nom, int nb_poi, int nb_caps, bool faisable,
                         int utilises, int nb_couverts, long long temps, int gen, const string& csv) {
    ofstream f(csv, ios::app);
    f << nom << "," << nb_poi << "," << nb_caps << "," << (faisable ? "Oui" : "Non") << ","
      << utilises << "," << nb_couverts << "," << temps << "," << gen << "\n"; // ajout nb_couverts
    f.close();
}


// Initialisation fichier CSV (sans écraser si déjà existant)
void initialiserFichierResultats(const string& nomFichier) {
    ifstream f_existant(nomFichier);
    if (!f_existant.good()) { 
        // Si le fichier n'existe pas, on le crée avec l'entête
        ofstream f(nomFichier);
        f << "Fichier,POI,Capteurs,Faisable,Capteurs_utilisés,POI_Couverts,Temps_ms,Generations\n";
        f.close();
    }
}
int compter_capteurs_utilises(const Chromosome& chromosome) {
    return count(chromosome.begin(), chromosome.end(), 1);
}

int compter_poi_couverts(const Chromosome& individu) {
    set<int> couverts;

    int nb_emplacements_par_capteur = N_EMPLACEMENTS;
    for (int capteur = 0; capteur < N_CAPTEURS; capteur++) {
        int debut = capteur * nb_emplacements_par_capteur;
        int id_emplacement = -1;

        for (int pos = debut; pos < debut + nb_emplacements_par_capteur; ++pos) {
            if (individu[pos] == 1) {
                id_emplacement = pos - debut;
                break;
            }
        }

        if (id_emplacement == -1) continue;

        for (size_t j = 0; j < points_interet.size(); ++j) {
            if (matrice_distance[id_emplacement][j] <= capteurs[capteur].rayon) {
                couverts.insert(j);
            }
        }
    }
    return couverts.size();
}


int main() {
    auto debut = chrono::high_resolution_clock::now();
    ofstream fichier("log.txt");
    fichier.close(); // vider log.txt

    initialiserFichierResultats("resultatsaz333.csv"); 
    lireJSON("config2.json");
    calculerMatriceDistance();
    calculerEtAfficherMatriceCouvertureParCapteur(); // afficher la matrice de chaque capteur

    auto population = initialiser_population();
    afficher_population(population);
    // Étapes de sélection, croisement et mutation AVANT la première génération
log("=== Sélection, Croisement et Mutation avant génération 0 ===");
cout << "=== Sélection, Croisement et Mutation avant génération 0 ===" << endl;

Population nouvelle_population;

while (nouvelle_population.size() < TAILLE_POPULATION) {
    Chromosome parent1 = tournoi_selection(population);
Chromosome parent2 = tournoi_selection(population);

    log("Parents sélectionnés :");
    string ligne1 = "P1: ";
    string ligne2 = "P2: ";
    for (int val : parent1) ligne1 += to_string(val) + " ";
    for (int val : parent2) ligne2 += to_string(val) + " ";
    log(ligne1);
    log(ligne2);

    Chromosome enfant1 = parent1;
    Chromosome enfant2 = parent2;

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0.0, 1.0);

 if (dis(gen) < CROSSOVER_RATE) {
    Chromosome e1, e2;
    int point;
    tie(e1, e2, point) = croisement(parent1, parent2);
    enfant1 = e1;
    enfant2 = e2;

    // Affichage plus complet :
    cout << "=== Croisement effectué ===" << endl;
    cout << "Point de croisement : " << point << endl;
    cout << "Parent 1: ";
    for (int val : parent1) cout << val << " ";
    cout << endl;
    cout << "Parent 2: ";
    for (int val : parent2) cout << val << " ";
    cout << endl;
    cout << "Enfant 1: ";
    for (int val : enfant1) cout << val << " ";
    cout << endl;
    cout << "Enfant 2: ";
    for (int val : enfant2) cout << val << " ";
    cout << endl;
    cout << "============================" << endl;
}


mutation(enfant1, nouvelle_population.size() + 1);
mutation(enfant2, nouvelle_population.size() + 2);


    log("Après mutation:");
    string mutation1 = "M1: ";
    string mutation2 = "M2: ";
    for (int val : enfant1) mutation1 += to_string(val) + " ";
    for (int val : enfant2) mutation2 += to_string(val) + " ";
    log(mutation1);
    log(mutation2);

    nouvelle_population.push_back(enfant1);
    if (nouvelle_population.size() < TAILLE_POPULATION)
        nouvelle_population.push_back(enfant2);
            cout << "============================" << endl;
}

log("=== Nouvelle population après sélection, croisement, mutation ===");
cout << "=== Nouvelle population après sélection, croisement, mutation ===" << endl;
afficher_population(nouvelle_population);

// Remplacer la population initiale par la nouvelle
population = nouvelle_population;


    Chromosome meilleur;
    double best_fitness = 0.0;
    int generation_finale = 0;
    bool faisable = false;




// Variables pour détecter la stagnation
int stagnation = 0;
int stagnation_max = 10; // Exemple: 10 générations sans amélioration
double best_fitness_precedent = 0.0;

for (int gen = 0; gen < N_GENERATIONS; ++gen) {
    log("*******************");
    log("Génération " + to_string(gen));
    log("*******************");

    for (const Chromosome& ind : population) {
        double fit = fonctionFitness(ind);

        string texte = "Individu : ";
        for (int gene : ind) texte += to_string(gene) + " ";
        texte += "| Fitness = " + to_string(fit);
        log(texte);

        if (fit > best_fitness) {
            best_fitness = fit;
            meilleur = ind;
            generation_finale = gen;
            faisable = true;
        }
    }

    log("Fin de l'évaluation de la génération.");
    log("--------------------------------------------------");

    // Détection stagnation
    if (best_fitness == best_fitness_precedent) {
        stagnation++;
    } else {
        stagnation = 0;
        best_fitness_precedent = best_fitness;
    }

    if (stagnation >= stagnation_max) {
        cout << "Arrêt anticipé : stagnation détectée (" << stagnation_max << " générations sans amélioration)." << endl;
        log("Arrêt anticipé : stagnation détectée.");
        break;
    }

    if (gen != N_GENERATIONS - 1) {
        population = nouvelle_generation(population);
    }
}


    if (faisable) {
        log("Meilleure solution trouvée à la génération " + to_string(generation_finale));
        log("Fitness : " + to_string(best_fitness));
        log("Chromosome : ");
        string texte = "";
        for (int gene : meilleur) texte += to_string(gene) + " ";
        log(texte);
    } else {
        log("Aucune solution réalisable trouvée.");
    }

    cout << "Fin de l'exécution. Voir log.txt." << endl;
     auto fin = chrono::high_resolution_clock::now();
     auto temps_ms = chrono::duration_cast<chrono::milliseconds>(fin - debut).count();
     cout << "Temps total d'exécution : " << temps_ms << " ms\n";

int meilleur_nb_capteurs = compter_capteurs_utilises(meilleur);

    cout << "\n=== Résultat final ===\n";
    if (faisable) {
        cout << "Solution trouvée en " << generation_finale << " générations.\n";
        cout << "Nombre de capteurs utilisés : " << meilleur_nb_capteurs << "\n";
        int nb_poi_couverts = compter_poi_couverts(meilleur);
        cout << "Points d'intérêt couverts : " << nb_poi_couverts << "\n";
        afficher_chromosome(meilleur);

        enregistrerResultat("config2.json", points_interet.size(), capteurs.size(), true,
                            meilleur_nb_capteurs, nb_poi_couverts, temps_ms, generation_finale, "resultatsaz333.csv");
    } else {
        cout << "Pas de solution faisable trouvée.\n";
        enregistrerResultat("config2.json", points_interet.size(), capteurs.size(), false,
                            0, 0, temps_ms, N_GENERATIONS, "resultatsaz333.csv");
    }
    return 0;
}
