#include <iostream> // Pour l'affichage console
#include <vector> // Pour utiliser les vecteurs dynamiques
#include <cmath> // Pour les fonctions mathématiques (pow, sqrt, etc.)
#include <fstream> // Pour la lecture et écriture dans des fichiers
#include <set>  // Pour gérer des ensembles sans doublons
#include <random> // Pour la génération de nombres aléatoires
//#include "json.hpp" // Librairie JSON pour lire la configuration depuis un fichier .json
#include <unordered_set>
#include <algorithm>
#include <iomanip> 
#include <nlohmann/json.hpp>
using namespace std;
using json = nlohmann::json; // Alias pour simplifier l'utilisation de la bibliothèque JSON


// === Structures de données utilisées ===
struct Position {
    double x;  // Coordonnée X dans le plan
    double y;  // Coordonnée Y dans le plan
};

struct PointInteret {
    int id;        // Identifiant unique
    Position pos;  // Position géographique du point
};

struct Capteur {
    int id;        // Identifiant du capteur
    double rayon;  // Rayon de couverture du capteur
};
struct Individu {
    std::vector<int> genes;// Positions des capteurs (-1 = désactivé)
    double fitness;//Score de qualité

    // Pour trier facilement les individus
    bool operator<(const Individu& other) const {
        return fitness > other.fitness; // Tri décroissant
    }
};

// === Variables globales de configuration ===
int TAILLE_POPULATION;       // Nombre d'individus dans chaque génération
int N_GENERATIONS;           // Nombre total de générations à simuler
double MUTATION_RATE;        // Taux de mutation appliqué aux chromosomes
double CROSSOVER_RATE;       // Taux de croisement entre individus
int TOURNAMENT_SIZE =0;         // Taille d’un tournoi lors de la sélection
int N_CAPTEURS = 0;          // Nombre de capteurs à placer
int N_EMPLACEMENTS = 0;      // Nombre de positions possibles pour chaque capteur
int chromosome_size = 0;     // Taille d'un chromosome (produit de capteurs × emplacements)
// === Conteneurs globaux ===
vector<Capteur> capteurs;                   // Liste de capteurs disponibles
vector<PointInteret> points_interet;        // Points d’intérêt à couvrir
vector<Position> emplacements;              // Positions où l’on peut placer un capteur
vector<vector<double>> matrice_distance;    // Matrice des distances entre emplacements et points

// === Alias pour types génétiques ===
typedef vector<int> Chromosome;             // 
typedef vector<Chromosome> Population;      // Ensemble d'individus formant une population
Chromosome meilleure_solution;
double meilleure_fitness = 0.0;

void reparer(Chromosome& chromo);

// Permet d'écrire une ligne de texte dans le fichier log.txt (mode append)
void log(string texte) {
    ofstream fichier("log.txt", ios::app); // ouverture en mode ajout
    if (fichier.is_open()) {
        fichier << texte << endl;          // écrit le texte dans le fichier
        fichier.close();                   // ferme le fichier après écriture
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
}// Chaque chromosome représente une configuration possible des capteurs

vector<vector<int>> initialiser_population() {
    vector<vector<int>> population;
    random_device rd;
    mt19937 gen(rd()); // Générateur aléatoire initialisé avec une graine système

    int nb_emplacements_par_capteur = N_EMPLACEMENTS;

    // Boucle de création des individus (chromosomes)
    for (int i = 0; i < TAILLE_POPULATION; i++) {
        vector<int> chromosome(N_CAPTEURS, -1); // Initialisé à -1 (capteur non placé)
        unordered_set<int> used; // Stocke les emplacements occupés (hors -1)

        for (int j = 0; j < N_CAPTEURS; j++) {
            uniform_int_distribution<> dis(-1, nb_emplacements_par_capteur - 1);
            int position;

            do {
                position = dis(gen); // Sélectionne un emplacement au hasard
            } while (position != -1 && used.count(position)); // Vérifie que l'emplacement est unique (sauf -1)

            chromosome[j] = position; // Assigne la position au capteur
            if (position != -1) {
                used.insert(position); // Ajoute l'emplacement aux utilisés
            }
        }

        population.push_back(chromosome); // Ajoute le chromosome à la population
    }

    return population;
}



double fonctionFitness(const Chromosome& individu);
// Affiche et loggue les individus de la population ainsi que leur fitness
void afficher_population(const vector<vector<int>>& population) {
    ofstream fichier("log.txt", ios::app); // ouvre le fichier log en mode ajout

    log("Population initiale :");
    for (size_t i = 0; i < population.size(); ++i) {
        string ligne = "Individu " + to_string(i + 1) + " : ";
      //  bool first = true;
        // Parcourt le chromosome pour afficher les positions activées 
        for (int j = 0; j < population[i].size(); j++) {
            ligne += to_string(population[i][j]);
            ligne += " , ";
           
        }
        // Calcule la fitness de l'individu
        double fitness = fonctionFitness(population[i]);
        ligne += " | Fitness = " + to_string(fitness);

        // Affiche sur la console
        cout << ligne << endl;

        // Écrit aussi dans le fichier log
        if (fichier.is_open()) {
            fichier << ligne << endl;
        }
    }

    log("--------------------------------------------------");
    fichier.close(); // ferme le fichier après écriture
}

// Vérifie si un individu couvre tous les points d'intérêt

//
bool estCouvertureComplete(const Chromosome& individu) {
    set<int> couverts;

    for (int capteur = 0; capteur < N_CAPTEURS; capteur++) {
        int id_emplacement = individu[capteur];
        if (id_emplacement == -1) continue;

        double rayon = capteurs[capteur].rayon;
        // Même si rayon == 0, le calcul se fera mais ne couvrira aucun point
        for (size_t j = 0; j < points_interet.size(); ++j) {
            if (matrice_distance[id_emplacement][j] <= rayon) {
                couverts.insert(j);
            }
        }
    }

    return couverts.size() == points_interet.size();
}

int CountCoveredPOI(const Chromosome& individu) {
    set<int> couverts;

    for (int capteur = 0; capteur < N_CAPTEURS; capteur++) {
        int id_emplacement = individu[capteur];
        if (id_emplacement == -1) continue;

        double rayon = capteurs[capteur].rayon;
        // Même si rayon == 0, le calcul se fera mais ne couvrira aucun point
        for (size_t j = 0; j < points_interet.size(); ++j) {
            if (matrice_distance[id_emplacement][j] <= rayon) {
                couverts.insert(j);
            }
        }
    }

    return couverts.size() ;
}

/*void eliminerCapteursInutiles(Chromosome& individu) {
    vector<bool> deja_couvert(points_interet.size(), false);

    // Marque les points déjà couverts par les capteurs utiles
    for (int capteur = 0; capteur < N_CAPTEURS; ++capteur) {
        int emplacement = individu[capteur];
        if (emplacement == -1) continue;

        double rayon = capteurs[capteur].rayon;
        for (size_t j = 0; j < points_interet.size(); ++j) {
            if (matrice_distance[emplacement][j] <= rayon) {
                deja_couvert[j] = true;
            }
        }
    }

    // Pour chaque capteur, on vérifie s’il est indispensable
    for (int capteur = 0; capteur < N_CAPTEURS; ++capteur) {
        int emplacement = individu[capteur];
        if (emplacement == -1) continue;

        double rayon = capteurs[capteur].rayon;
        bool est_util = false;

        // Simuler la suppression du capteur
        for (size_t j = 0; j < points_interet.size(); ++j) {
            if (matrice_distance[emplacement][j] <= rayon) {
                // On vérifie si ce point serait découvert sans ce capteur
                bool couvert_autrement = false;
                for (int autre = 0; autre < N_CAPTEURS; ++autre) {
                    if (autre == capteur) continue;
                    int autre_emp = individu[autre];
                    if (autre_emp == -1) continue;

                    double r2 = capteurs[autre].rayon;
                    if (matrice_distance[autre_emp][j] <= r2) {
                        couvert_autrement = true;
                        break;
                    }
                }
                if (!couvert_autrement) {
                    est_util = true;
                    break; // Dès qu’un point n’est couvert que par ce capteur
                }
            }
        }

        if (!est_util) {
            individu[capteur] = -1; // Marquer ce capteur comme inutile
        }
    }
}*/
void reparer(Chromosome& chromo) {
    set<int> emplacements_utilises;
    vector<bool> points_couverts(points_interet.size(), false);
    int capteurs_utilises = 0;

    // Étape 1: Réparer les doublons intelligemment
    for (int capteur = 0; capteur < N_CAPTEURS; ++capteur) {
        if (chromo[capteur] != -1) {
            if (emplacements_utilises.count(chromo[capteur])) {
                // Trouver le meilleur emplacement alternatif
                int meilleur_emp = -1;
                double meilleure_couverture = -1;
                
                for (int emp = 0; emp < N_EMPLACEMENTS; ++emp) {
                    if (!emplacements_utilises.count(emp)) {
                        double couverture = 0;
                        for (size_t j = 0; j < points_interet.size(); ++j) {
                            if (!points_couverts[j] && matrice_distance[emp][j] <= capteurs[capteur].rayon) {
                                couverture++;
                            }
                        }
                        if (couverture > meilleure_couverture) {
                            meilleure_couverture = couverture;
                            meilleur_emp = emp;
                        }
                    }
                }
                
                if (meilleur_emp != -1) {
                    chromo[capteur] = meilleur_emp;
                } else {
                    chromo[capteur] = -1; // Désactiver ce capteur
                    continue;
                }
            }
            
            // Mettre à jour la couverture
            int emp = chromo[capteur];
            emplacements_utilises.insert(emp);
            capteurs_utilises++;
            double rayon = capteurs[capteur].rayon;
            for (size_t j = 0; j < points_interet.size(); ++j) {
                if (matrice_distance[emp][j] <= rayon) {
                    points_couverts[j] = true;
                }
            }
        }
    }

    // Étape 2: Amélioration progressive plutôt que couverture totale obligatoire
    bool amelioration;
    do {
        amelioration = false;
        
        // Essayer d'ajouter/améliorer un capteur
        for (int capteur = 0; capteur < N_CAPTEURS; ++capteur) {
            if (chromo[capteur] == -1) { // Capteur disponible
                // Trouver le meilleur emplacement pour ce capteur
                int meilleur_emp = -1;
                int max_couverture = 0;
                
                for (int emp = 0; emp < N_EMPLACEMENTS; ++emp) {
                    if (!emplacements_utilises.count(emp)) {
                        int couverture = 0;
                        for (size_t j = 0; j < points_interet.size(); ++j) {
                            if (!points_couverts[j] && matrice_distance[emp][j] <= capteurs[capteur].rayon) {
                                couverture++;
                            }
                        }
                        if (couverture > max_couverture) {
                            max_couverture = couverture;
                            meilleur_emp = emp;
                        }
                    }
                }
                
                if (meilleur_emp != -1 && max_couverture > 0) {
                    chromo[capteur] = meilleur_emp;
                    emplacements_utilises.insert(meilleur_emp);
                    capteurs_utilises++;
                    double rayon = capteurs[capteur].rayon;
                    for (size_t j = 0; j < points_interet.size(); ++j) {
                        if (matrice_distance[meilleur_emp][j] <= rayon) {
                            points_couverts[j] = true;
                        }
                    }
                    amelioration = true;
                }
            }
        }
    } while (amelioration);
}

// Évalue la qualité (fitness) d'un individu (chromosome) en fonction de la couverture et du nombre de capteurs utiles
double fonctionFitness(const Chromosome& individu) {
    
    int capteurs_utilises = count_if(individu.begin(), individu.end(), 
                                   [](int emp) { return emp != -1; });
    
    // Priorité 1: couverture complète
    // Priorité 2: minimiser le nombre de capteurs
    return pow(CountCoveredPOI(individu),2) / pow((1.0 + capteurs_utilises),4); 
}
// Doit être déclarée ailleurs
bool estCouvertureComplete(const vector<int>& individu);
double fonctionFitness(const vector<int>& individu);

void sauvegarderMeilleurIndividu(const vector<vector<int>>& population, const string& nomFichier) {
    double meilleureFitness = -1.0;
    vector<int> meilleurIndividu;

    for (const auto& individu : population) {
        if (estCouvertureComplete(individu)) {
            double fitness = fonctionFitness(individu);
            if (meilleureFitness < 0 || fitness > meilleureFitness) {
                meilleureFitness = fitness;
                meilleurIndividu = individu;
            }
        }
    }

    ofstream fichier(nomFichier);
    if (fichier.is_open()) {
        if (meilleurIndividu.empty()) {
            fichier << "Aucune solution valide trouvée dans cette génération.\n";
            // Ajout pour affichage console
            cout << "Aucune solution valide dans la génération actuelle.\n";
        } else {
            fichier << "Meilleure solution (fitness = " << meilleureFitness << ") :\n";
            fichier << "Nombre de capteurs utilisés: " 
                   << count_if(meilleurIndividu.begin(), meilleurIndividu.end(), 
                             [](int x) { return x != -1; }) << "\n";
            fichier << "Configuration: ";
            for (int v : meilleurIndividu) {
                fichier << v << " ";
            }
            fichier << endl;
            
            // Affichage console
            cout << "Meilleure solution de la génération (fitness = " 
                 << meilleureFitness << "):\n";
            for (size_t i = 0; i < meilleurIndividu.size(); ++i) {
                if (meilleurIndividu[i] != -1) {
                    cout << "Capteur " << i << " -> Position " 
                         << meilleurIndividu[i] << endl;
                }
            }
        }
        fichier.close();
    }
}
vector<vector<int>> selectionParTournoiAvecLog(const vector<vector<int>>& population, int tournament_size) {
    ofstream fichier("log.txt", ios::app); // ouvre le fichier log en mode ajout
    if (!fichier.is_open()) {
        throw runtime_error("Erreur lors de l'ouverture de log.txt");
    }

    fichier << "=== Début de la sélection par tournoi ===\n";

    vector<vector<int>> nouvellePopulation;
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dist(0, population.size() - 1);

    for (size_t i = 0; i < population.size(); ++i) {
        double meilleureFitness = -1.0;
        vector<int> meilleurIndividu;

        fichier << "Tournoi #" << i + 1 << ":\n";

        for (int j = 0; j < tournament_size; ++j) {
            const vector<int>& candidat = population[dist(gen)];
            double fitness = fonctionFitness(candidat);

            fichier << "  Candidat " << j + 1 << " : [";
            for (int g : candidat) fichier << g << " ";
            fichier << "] fitness = " << fitness << "\n";

            if (meilleurIndividu.empty() || fitness > meilleureFitness) {
                meilleureFitness = fitness;
                meilleurIndividu = candidat;
            }
        }

        fichier << "  → Sélectionné : [";
        for (int g : meilleurIndividu) fichier << g << " ";
        fichier << "] fitness = " << meilleureFitness << "\n\n";

        nouvellePopulation.push_back(meilleurIndividu);
    }

    fichier << "=== Fin de la sélection ===\n";
    fichier << "--------------------------------------------------\n";
    fichier.close(); // ferme correctement le fichier

    return nouvellePopulation;
}



// Lit la configuration du problème depuis un fichier JSON et initialise les structures globales
void lireJSON(string nom_fichier) {
    ifstream input(nom_fichier); // Ouvre le fichier JSON
    json j;
    input >> j; // Parse le fichier dans un objet JSON

    // Chargement des capteurs
    for (auto& c : j["capteurs"]) {
        capteurs.push_back({c["id"], c["rayon"]}); 
    }
    // Chargement des points d'intérêt
    for (auto& p : j["points_interet"]) {
        points_interet.push_back({p["id"], {p["x"], p["y"]}});
    }

    // Chargement des emplacements
    for (auto& e : j["emplacements"]) {
        emplacements.push_back({e["x"], e["y"]});
    }

    // Chargement des paramètres génétiques
    TAILLE_POPULATION = j["genetic_params"]["population_size"];
    N_GENERATIONS = j["genetic_params"]["generations"];
    MUTATION_RATE = j["genetic_params"]["mutation_rate"];
    CROSSOVER_RATE = j["genetic_params"]["crossover_rate"];
    TOURNAMENT_SIZE = j["genetic_params"]["tournament_size"];

    // Calculs dérivés
    N_CAPTEURS = capteurs.size();
    N_EMPLACEMENTS = emplacements.size();
    chromosome_size = N_CAPTEURS ; // Longueur totale du chromosome
} 
//
// Corriger les doublons dans un chromosome
Chromosome fixDuplicates(const Chromosome& chromosome, int chromosome_size) {
    set<int> seen;
    set<int> all_genes;
    for (int i = 0; i < chromosome_size; ++i) {
        all_genes.insert(i);
    }

    Chromosome fixed;
    for (int gene : chromosome) {
        if (seen.find(gene) == seen.end()) {
            fixed.push_back(gene);
            seen.insert(gene);
        } else {
            // Trouver un gène non utilisé
            set<int> unused_genes;
            set_difference(all_genes.begin(), all_genes.end(),
                           seen.begin(), seen.end(),
                           inserter(unused_genes, unused_genes.begin()));

            int replacement = *unused_genes.begin();
            fixed.push_back(replacement);
            seen.insert(replacement);
        }
    }
    return fixed;
}

// Fonction de croisement
Population crossover(const Population& selected_population, float crossover_rate, int chromosome_size) {
    ofstream fichier("log.txt", ios::app);
    Population new_population;
    vector<Chromosome> parents = selected_population;

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> prob(0.0, 1.0);
    shuffle(parents.begin(), parents.end(), gen);

    fichier << "\n=== DEBUT CROISEMENT ===\n";
    fichier << "Nombre de parents: " << parents.size() << "\n";
    fichier << "Taux de croisement: " << crossover_rate << "\n";

    for (size_t i = 0; i + 1 < parents.size(); i += 2) {
        fichier << "\nPaire de parents " << i/2 + 1 << ":\n";
        fichier << "Parent 1: [ ";
        for (int g : parents[i]) fichier << g << " ";
        fichier << "] Fitness: " << fonctionFitness(parents[i]) << "\n";
        fichier << "Parent 2: [ ";
        for (int g : parents[i+1]) fichier << g << " ";
        fichier << "] Fitness: " << fonctionFitness(parents[i+1]) << "\n";

        if (prob(gen) < crossover_rate) {
            uniform_int_distribution<> point_dist(1, chromosome_size - 1);
            int crossover_point = point_dist(gen);
            
            Chromosome child1, child2;
            
            // Premier enfant
            child1.insert(child1.end(), parents[i].begin(), parents[i].begin() + crossover_point);
            child1.insert(child1.end(), parents[i+1].begin() + crossover_point, parents[i+1].end());
            
            // Deuxième enfant
            child2.insert(child2.end(), parents[i+1].begin(), parents[i+1].begin() + crossover_point);
            child2.insert(child2.end(), parents[i].begin() + crossover_point, parents[i].end());
            
            fichier << "Croisement au point " << crossover_point << "\n";
            fichier << "Enfant 1 avant réparation: [ ";
            for (int g : child1) fichier << g << " ";
            fichier << "]\n";
            fichier << "Enfant 2 avant réparation: [ ";
            for (int g : child2) fichier << g << " ";
            fichier << "]\n";

            reparer(child1);
            reparer(child2);
            
            fichier << "Enfant 1 après réparation: [ ";
            for (int g : child1) fichier << g << " ";
            fichier << "] Fitness: " << fonctionFitness(child1) << "\n";
            fichier << "Enfant 2 après réparation: [ ";
            for (int g : child2) fichier << g << " ";
            fichier << "] Fitness: " << fonctionFitness(child2) << "\n";
            
            new_population.push_back(child1);
            new_population.push_back(child2);
        } else {
            fichier << "Pas de croisement (taux non atteint)\n";
            fichier << "Enfants identiques aux parents\n";
            new_population.push_back(parents[i]);
            new_population.push_back(parents[i+1]);
        }
    }
    
    if (parents.size() % 2 != 0) {
        fichier << "\nParent impair non apparié ajouté tel quel: [ ";
        for (int g : parents.back()) fichier << g << " ";
        fichier << "] Fitness: " << fonctionFitness(parents.back()) << "\n";
        new_population.push_back(parents.back());
    }
    
    fichier << "=== FIN CROISEMENT ===\n";
    fichier << "Nombre total d'enfants générés: " << new_population.size() << "\n\n";
    fichier.close();
    
    return new_population;
}
// Applique une mutation sur chaque individu de la population
// Applique une mutation sur chaque individu de la population avec log détaillé
void mutation(Population& population, double mutation_rate, int n_emplacements) {
    ofstream fichier("log.txt", ios::app);
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0.0, 1.0);
    
    fichier << "\n=== DEBUT MUTATION ===\n";
    fichier << "Taux de mutation: " << mutation_rate << "\n";
    fichier << "Nombre d'individus: " << population.size() << "\n";
    
    int total_mutations = 0;
    
    for (size_t idx = 0; idx < population.size(); ++idx) {
        bool individu_mute = false;
        fichier << "\nIndividu " << idx+1 << " avant mutation: [ ";
        for (int g : population[idx]) fichier << g << " ";
        fichier << "] Fitness: " << fonctionFitness(population[idx]) << "\n";
        
        for (auto& gene : population[idx]) {
            if (dis(gen) < mutation_rate) {
                int old_gene = gene;
                if (dis(gen) < 0.5) {
                    gene = -1; //-1 semble représenter un capteur désactivé/Cela permet d'explorer des solutions avec moins de capteurs
                } else {
                    gene = rand() % n_emplacements;
                }
                fichier << "Mutation: " << old_gene << " -> " << gene << "\n";
                total_mutations++;
                individu_mute = true;
            }
        }
        
        if (individu_mute) {
            fichier << "Individu " << idx+1 << " après mutation: [ ";
            for (int g : population[idx]) fichier << g << " ";
            fichier << "]\n";
            
            reparer(population[idx]);
            fichier << "Après réparation: [ ";
            for (int g : population[idx]) fichier << g << " ";
            fichier << "] Nouvelle fitness: " << fonctionFitness(population[idx]) << "\n";
        } else {
            fichier << "Aucune mutation pour cet individu\n";
        }
    }
    
    fichier << "\nTotal mutations appliquées: " << total_mutations << "\n";
    fichier << "=== FIN MUTATION ===\n\n";
    fichier.close();
}
void remplacement(Population& population, const Population& enfants, int elite_count) {
    ofstream fichier("log.txt", ios::app);
    
    fichier << "\n=== DEBUT REMPLACEMENT ===\n";
    fichier << "Taille population actuelle: " << population.size() << "\n";
    fichier << "Nombre d'enfants: " << enfants.size() << "\n";
    fichier << "Nombre d'élites à conserver: " << elite_count << "\n";
    
    // Trier la population actuelle
    sort(population.begin(), population.end(), [](const Chromosome& a, const Chromosome& b) {
        return fonctionFitness(a) > fonctionFitness(b);
    });
    
    fichier << "\nMeilleurs individus actuels:\n";
    for (int i = 0; i < min(3, (int)population.size()); ++i) {
        fichier << i+1 << ". [ ";
        for (int g : population[i]) fichier << g << " ";
        fichier << "] Fitness: " << fonctionFitness(population[i]) << "\n";
    }
    
    // Trier les enfants
    Population enfants_tries = enfants;
    sort(enfants_tries.begin(), enfants_tries.end(), [](const Chromosome& a, const Chromosome& b) {
        return fonctionFitness(a) > fonctionFitness(b);
    });
    
    fichier << "\nMeilleurs enfants:\n";
    for (int i = 0; i < min(3, (int)enfants_tries.size()); ++i) {
        fichier << i+1 << ". [ ";
        for (int g : enfants_tries[i]) fichier << g << " ";
        fichier << "] Fitness: " << fonctionFitness(enfants_tries[i]) << "\n";
    }
    
    // Nouvelle population
    Population nouvelle_population;
    
    // Élitisme
    for (int i = 0; i < elite_count && i < population.size(); ++i) {
        nouvelle_population.push_back(population[i]);
    }
    fichier << "\nÉlites conservées: " << nouvelle_population.size() << "\n";
    
    // Ajout des meilleurs enfants
    int index_enfant = 0;
    while (nouvelle_population.size() < population.size() &&
       static_cast<size_t>(index_enfant) < enfants_tries.size())
{
        nouvelle_population.push_back(enfants_tries[index_enfant]);
        index_enfant++;
    }
    fichier << "Enfants ajoutés: " << index_enfant << "\n";
    
    // Complément avec nouveaux individus si nécessaire
    if (nouvelle_population.size() < population.size()) {
        int manquants = population.size() - nouvelle_population.size();
        fichier << "Génération de " << manquants << " nouveaux individus\n";
        
        for (int i = 0; i < manquants; ++i) {
            Chromosome nouveau;
            for (int j = 0; j < N_CAPTEURS; ++j) {
                nouveau.push_back(rand() % N_EMPLACEMENTS);
            }
            reparer(nouveau);
            nouvelle_population.push_back(nouveau);
            
            fichier << "Nouvel individu " << i+1 << ": [ ";
            for (int g : nouveau) fichier << g << " ";
            fichier << "] Fitness: " << fonctionFitness(nouveau) << "\n";
        }
    }
    
    population = nouvelle_population;
    
    fichier << "\nNouvelle population:\n";
    double avg_fitness = 0;
    double max_fitness = 0;
    for (const auto& indiv : population) {
        double fit = fonctionFitness(indiv);
        avg_fitness += fit;
        if (fit > max_fitness) max_fitness = fit;
    }
    avg_fitness /= population.size();
    
    fichier << "Fitness max: " << max_fitness << "\n";
    fichier << "Fitness moyenne: " << avg_fitness << "\n";
    fichier << "=== FIN REMPLACEMENT ===\n\n";
    fichier.close();
}
/*

}*/
int main() {
    // Vider le fichier log
    ofstream fichier("log.txt");
    fichier.close();

    // Charger les données
    lireJSON("config1.json");
    calculerMatriceDistance();
    calculerEtAfficherMatriceCouvertureParCapteur();

    // Initialisation population
    Population population = initialiser_population();
    afficher_population(population);

    int generation = 0;
    int stagnation_count = 0;
    const int MAX_GENERATIONS = 10;  // Augmenté
    const int MAX_STAGNATION = 10;    // Augmenté

    double meilleure_fitness_precedente = -1.0;
    Chromosome meilleure_solution_globale;
    double meilleure_fitness_globale = -1.0;

    // === BOUCLE PRINCIPALE ===
    while (generation <  MAX_GENERATIONS && stagnation_count <  MAX_STAGNATION) {
        log("\n--- Génération " + to_string(generation + 1) + " ---");

         // Étape 1 : Évaluation
    int solutions_valides = 0; // Compte les individus qui couvrent tous les points
    vector<double> fitness_values; // Stocke les valeurs de fitness de chaque individu
    for (const auto& individu : population) {    // Parcours de chaque individu de la population actuelle
        double fitness = fonctionFitness(individu); // Calcul de la fitness pour l'individu courant
        fitness_values.push_back(fitness);// Enregistrement de cette fitness dans le tableau des valeurs
        if (fitness > 0) solutions_valides++; // Si la fitness > 0, alors l'individu couvre tous les points (c'est une solution valide)
        
        // Mise à jour meilleure solution globale...
    }
    
    // Affichage du nombre de solutions valides  // Affichage console (et éventuellement log) du nombre de solutions valides pour cette génération
    cout << "Génération " << generation + 1 
         << " - Solutions valides: " << solutions_valides 
         << "/" << population.size() << endl;
        // Étape 1 : Évaluation
      //  vector<double> fitness_values;
        for (const auto& individu : population) {// Ici, on cherche le meilleur individu global (parmi toutes les générations passées)
            double fitness = fonctionFitness(individu);
            fitness_values.push_back(fitness);
            
            // Mise à jour meilleure solution globale// et que l'individu couvre tous les points (vérification de validité)
        // Si on trouve une fitness supérieure à celle enregistrée jusque-là,
            if (fitness > meilleure_fitness_globale && estCouvertureComplete(individu)) {
                meilleure_fitness_globale = fitness;  // On met à jour la meilleure fitness connue globalement
                meilleure_solution_globale = individu;  // On sauvegarde l'individu correspondant comme meilleure solution
                //eliminerCapteursInutiles(meilleure_solution_globale); // Optimisation immédiate // On élimine les capteurs inutiles pour simplifier la solution
            }
        }

        // Étape 2 : Sélection 
        Population selected_population = selectionParTournoiAvecLog(population, TOURNAMENT_SIZE);

        // Étape 3 : Croisement
        Population crossed_population = crossover(selected_population, CROSSOVER_RATE, chromosome_size);

        // Étape 4 : Mutation
        mutation(crossed_population, MUTATION_RATE, N_EMPLACEMENTS);

        // Étape 5 : Réparation des solutions invalides
        for (auto& individu : crossed_population) {
            if (!estCouvertureComplete(individu)) {
                reparer(individu);
            }
        }

        // Étape 6 : Remplacement avec élitisme
        remplacement(population, crossed_population, 2); // On garde les 2 meilleurs

        // Détection de stagnation
        double max_fitness = *max_element(fitness_values.begin(), fitness_values.end());
        if (max_fitness <= meilleure_fitness_precedente) {
            stagnation_count++;
            log("Stagnation #" + to_string(stagnation_count));
        } else {
            stagnation_count = 0;
            meilleure_fitness_precedente = max_fitness;
        }
        generation++;//Cela prépare la boucle à traiter la génération suivante.
    }

    
    // === FIN ALGORITHME ===
log("\n=== RESULTAT FINAL ===");
log("Nombre total de générations: " + to_string(generation));
log("Dernier compteur de stagnation: " + to_string(stagnation_count));

// Sauvegarde de la meilleure solution trouvée
ofstream final_solution("meilleure_solution_finale.txt");
if (final_solution.is_open()) {
    if (meilleure_fitness_globale <= 0) { // Aucune solution valide trouvée
        final_solution << "Aucune solution valide n'a été trouvée après " 
                      << generation << " générations.\n";
        cout << "\n=== AUCUNE SOLUTION VALIDE TROUVÉE ===" << endl;
        cout << "L'algorithme n'a pas trouvé de configuration de capteurs" << endl;
        cout << "couvrant tous les points d'intérêt après " << generation 
             << " générations." << endl;
    } else {
        final_solution << "Meilleure solution globale (fitness = " 
                      << meilleure_fitness_globale << "):\n";
        final_solution << "Nombre de capteurs utilisés: " 
                      << count_if(meilleure_solution_globale.begin(), 
                                 meilleure_solution_globale.end(), 
                                 [](int x) { return x != -1; }) << "\n";
        
        final_solution << "Configuration: ";
        for (int gene : meilleure_solution_globale) {
            final_solution << gene << " ";
        }
        
        // Affichage console
        cout << "\n=== SOLUTION OPTIMALE TROUVEE ===" << endl;
        cout << "Capteurs utilises: ";
        for (size_t i = 0; i < meilleure_solution_globale.size(); ++i) {
            if (meilleure_solution_globale[i] != -1) {
                cout << "Capteur " << i << " -> Position " 
                     << meilleure_solution_globale[i] << endl;
            }
        }
        cout << "Fitness: " << meilleure_fitness_globale << endl;
    }
    final_solution.close();
}
    // Affichage console
    cout << "\n=== SOLUTION OPTIMALE TROUVEE ===" << endl;
    cout << "Capteurs utilises: ";
    for (size_t i = 0; i < meilleure_solution_globale.size(); ++i) {
        if (meilleure_solution_globale[i] != -1) {
            cout << "Capteur " << i << " -> Position " 
                 << meilleure_solution_globale[i] << endl;
        }
    }
    cout << "Fitness: " << meilleure_fitness_globale << endl;

    return 0;
}
