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
#include <map>

#include <chrono>  // pour mesurer le temps

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
int TOURNAMENT_SIZE =0;         // Taille d'un tournoi lors de la sélection
int N_CAPTEURS = 0;          // Nombre de capteurs à placer
int N_EMPLACEMENTS = 0;      // Nombre de positions possibles pour chaque capteur
int chromosome_size = 0;     // Taille d'un chromosome (produit de capteurs × emplacements)
// === Conteneurs globaux ===
vector<Capteur> capteurs;                   // Liste de capteurs disponibles
vector<PointInteret> points_interet;        // Points d'intérêt à couvrir
vector<Position> emplacements;              // Positions où l'on peut placer un capteur
vector<vector<double>> matrice_distance;    // Matrice des distances entre emplacements et points

// === Alias pour types génétiques ===
typedef vector<int> Chromosome;             // 
typedef vector<Chromosome> Population;      // Ensemble d'individus formant une population
Chromosome meilleure_solution;
double meilleure_fitness = 0.0;

// Déclarations des fonctions
int CountCoveredPOI(const Chromosome& individu);
bool estCouvertureComplete(const Chromosome& individu);
double fonctionFitness(const Chromosome& individu);
void lireJSON(string nom_fichier, int instance_index);

//void reparer(Chromosome& chromo);

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

// Évalue la qualité (fitness) d'un individu (chromosome) en fonction de la couverture et du nombre de capteurs utiles
double fonctionFitness(const Chromosome& individu) {
    
    int capteurs_utilises = count_if(individu.begin(), individu.end(), 
                                   [](int emp) { return emp != -1; });//combien de capteur utilser
    
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
    // Initialise la nouvelle population qui sera retournée
    vector<vector<int>> nouvellePopulation;
    random_device rd;    // Initialise le générateur de nombres aléatoires
    mt19937 gen(rd());
    uniform_int_distribution<> dist(0, population.size() - 1);    // Distribution uniforme pour sélectionner aléatoirement des individus
    // Pour chaque individu à sélectionner dans la nouvelle population
    for (size_t i = 0; i < population.size(); ++i) {
        double meilleureFitness = -1.0;
        vector<int> meilleurIndividu;
        // Loggue le numéro du tournoi
        fichier << "Tournoi #" << i + 1 << ":\n";
        // Effectue le tournoi avec 'tournament_size' participants
        for (int j = 0; j < tournament_size; ++j) {
            const vector<int>& candidat = population[dist(gen)];// sélectionne un candidat aléatoire
            double fitness = fonctionFitness(candidat);// Calcule sa fitness
            // Loggue les détails du candidat
            fichier << "  Candidat " << j + 1 << " : [";
            for (int g : candidat) fichier << g << " ";
            fichier << "] fitness = " << fitness << "\n";
            // Met à jour le meilleur candidat du tournoi
            if (meilleurIndividu.empty() || fitness > meilleureFitness) {
                meilleureFitness = fitness;
                meilleurIndividu = candidat;
            }
        }
        // Loggue le gagnant du tournoi
        fichier << "  → Sélectionné : [";
        for (int g : meilleurIndividu) fichier << g << " ";
        fichier << "] fitness = " << meilleureFitness << "\n\n";
        // Ajoute le gagnant à la nouvelle population
        nouvellePopulation.push_back(meilleurIndividu);
    }
    // Loggue la fin de la sélection
    fichier << "=== Fin de la sélection ===\n";
    fichier << "--------------------------------------------------\n";
    fichier.close(); // ferme correctement le fichier
    afficher_population(nouvellePopulation);    // Affiche la nouvelle population

    return nouvellePopulation;
}



// Nouvelle fonction pour lire une instance depuis un objet json image
void lireInstanceDepuisObjet(const json& image_obj, const json& genetic_params) {
    // Chargement des capteurs
    for (auto& c : image_obj["sensors"]) {
        capteurs.push_back({c["id"], c["range"]});
    }
    // Chargement des points d'intérêt (POI)
    for (auto& p : image_obj["pois"]) {
        points_interet.push_back({0, {p["x"], p["y"]}});
    }
    // Chargement des emplacements
    for (auto& p : image_obj["locations"]) {
        emplacements.push_back({p["x"], p["y"]});
    }
    // Chargement des paramètres génétiques
    TAILLE_POPULATION = genetic_params["population_size"];
    N_GENERATIONS = genetic_params["generations"];
    MUTATION_RATE = genetic_params["mutation_rate"];
    CROSSOVER_RATE = genetic_params["crossover_rate"];
    TOURNAMENT_SIZE = genetic_params["tournament_size"];
    N_CAPTEURS = capteurs.size();
    N_EMPLACEMENTS = emplacements.size();
    chromosome_size = N_CAPTEURS;
}

// Corriger les doublons dans un chromosome
Chromosome fixDuplicates(const Chromosome& chromosome, int chromosome_size) {
    set<int> seen; // Ensemble pour suivre les gènes déjà vus
    set<int> all_genes;// Ensemble de tous les gènes possibles
    for (int i = 0; i < chromosome_size; ++i) {
        all_genes.insert(i);
    }
    Chromosome fixed;
        // Pour chaque gène dans le chromosome
    for (int gene : chromosome) {
        if (seen.find(gene) == seen.end()) { // Si le gène n'a pas encore été vu
            fixed.push_back(gene);
            seen.insert(gene);
        } else {
            // Trouver un gène non utilisé
            set<int> unused_genes;
            set_difference(all_genes.begin(), all_genes.end(),
                           seen.begin(), seen.end(),
                           inserter(unused_genes, unused_genes.begin()));
            // Prend le premier gène non utilisé
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
    vector<Chromosome> parents = selected_population;  // Copie la population sélectionnée
    // Initialise le générateur aléatoire
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> prob(0.0, 1.0);  // Pour décider si on fait un croisement
    //Elle mélange (shuffle) les chromosomes de la population sélectionnée de manière aléatoire
    // avant de former les paires de parents.
    shuffle(parents.begin(), parents.end(), gen);// Mélange les parents aléatoirement

    fichier << "\n=== DEBUT CROISEMENT ===\n";
    fichier << "Nombre de parents: " << parents.size() << "\n";
    fichier << "Taux de croisement: " << crossover_rate << "\n";
    // Parcourt les parents par paires
    for (size_t i = 0; i + 1 < parents.size(); i += 2) {
        fichier << "\nPaire de parents " << i/2 + 1 << ":\n";
        fichier << "Parent 1: [ ";        // Loggue le premier parent
        for (int g : parents[i]) fichier << g << " ";
        fichier << "] Fitness: " << fonctionFitness(parents[i]) << "\n";
        fichier << "Parent 2: [ ";        // Loggue le second parent
        for (int g : parents[i+1]) fichier << g << " ";
        fichier << "] Fitness: " << fonctionFitness(parents[i+1]) << "\n";
        // Décide si on fait un croisement
        if (prob(gen) < crossover_rate) {
                        // Choisit un point de croisement aléatoire
            uniform_int_distribution<> point_dist(1, chromosome_size - 1);
            int crossover_point = point_dist(gen); 
            Chromosome child1, child2;
            
        // Premier enfant: début parent1 + fin parent2
            child1.insert(child1.end(), parents[i].begin(), parents[i].begin() + crossover_point);
            child1.insert(child1.end(), parents[i+1].begin() + crossover_point, parents[i+1].end());
     // Second enfant: début parent2 + fin parent1
            child2.insert(child2.end(), parents[i+1].begin(), parents[i+1].begin() + crossover_point);
            child2.insert(child2.end(), parents[i].begin() + crossover_point, parents[i].end());
    // Loggue les enfants avant réparation  
            fichier << "Croisement au point " << crossover_point << "\n";
            fichier << "Enfant 1 avant réparation: [ ";
            for (int g : child1) fichier << g << " ";
            fichier << "]\n";
            fichier << "Enfant 2 avant réparation: [ ";
            for (int g : child2) fichier << g << " ";
            fichier << "]\n";
            // Répare les doublons
            child1 = fixDuplicates(child1, chromosome_size);
            child2 = fixDuplicates(child2, chromosome_size);
         // Loggue les enfants après réparation
            fichier << "Enfant 1 après réparation: [ ";
            for (int g : child1) fichier << g << " ";
            fichier << "] Fitness: " << fonctionFitness(child1) << "\n";
            fichier << "Enfant 2 après réparation: [ ";
            for (int g : child2) fichier << g << " ";
            fichier << "] Fitness: " << fonctionFitness(child2) << "\n";
     // Ajoute les enfants à la nouvelle population    
            new_population.push_back(child1);
            new_population.push_back(child2);
        } else {  // Si pas de croisement, ajoute les parents tels quels
            fichier << "Pas de croisement (taux non atteint)\n";
            fichier << "Enfants identiques aux parents\n";
            new_population.push_back(parents[i]);
            new_population.push_back(parents[i+1]);
        }
    }
// Gère le cas où le nombre de parents est impair
    if (parents.size() % 2 != 0) {
        fichier << "\nParent impair non apparié ajouté tel quel: [ ";
        for (int g : parents.back()) fichier << g << " ";
        fichier << "] Fitness: " << fonctionFitness(parents.back()) << "\n";
        new_population.push_back(parents.back());
    }
       // Loggue la fin du croisement 
    fichier << "=== FIN CROISEMENT ===\n";
    fichier << "Nombre total d'enfants générés: " << new_population.size() << "\n\n";
    fichier.close();
    afficher_population(new_population);   // Affiche la nouvelle population

    return new_population;
}
// Applique une mutation sur chaque individu de la population
// Applique une mutation sur chaque individu de la population avec log détaillé
void mutation(Population& population, double mutation_rate, int n_emplacements) {
    ofstream fichier("log.txt", ios::app);
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0.0, 1.0); // Pour décider si on mute un gène
       // Loggue le début de la mutation 
    fichier << "\n=== DEBUT MUTATION ===\n";
    fichier << "Taux de mutation: " << mutation_rate << "\n";
    fichier << "Nombre d'individus: " << population.size() << "\n";
    
    int total_mutations = 0;
       // Pour chaque individu de la population 
    for (size_t idx = 0; idx < population.size(); ++idx) {
        bool individu_mute = false;
        fichier << "\nIndividu " << idx+1 << " avant mutation: [ ";        // Loggue l'individu avant mutation
        for (int g : population[idx]) fichier << g << " ";
        fichier << "] Fitness: " << fonctionFitness(population[idx]) << "\n";
         // Pour chaque gène de l'individu
        for (auto& gene : population[idx]) {
            if (dis(gen) < mutation_rate) {  // Décide si on mute ce gène
                int old_gene = gene;
                if (dis(gen) < 0.5) { // Décide aléatoirement le type de mutation
                    gene = -1; //-1 semble représenter un capteur désactivé/Cela permet d'explorer des solutions avec moins de capteurs
                } else {
                    gene = rand() % n_emplacements;// Change la position
                }
                fichier << "Mutation: " << old_gene << " -> " << gene << "\n"; // Loggue la mutation
                total_mutations++;
                individu_mute = true;
            }

        }
                // Loggue l'individu après mutation si nécessaire
        if (individu_mute) {
            fichier << "Individu " << idx+1 << " après mutation: [ ";
            for (int g : population[idx]) fichier << g << " ";
            fichier << "]\n";
            
        //    reparer(population[idx]);
            fichier << "Après réparation: [ ";
            for (int g : population[idx]) fichier << g << " ";
            fichier << "] Nouvelle fitness: " << fonctionFitness(population[idx]) << "\n";
        } else {
            fichier << "Aucune mutation pour cet individu\n";
        }
    }
        // Loggue le nombre total de mutations
    fichier << "\nTotal mutations appliquées: " << total_mutations << "\n";
    fichier << "=== FIN MUTATION ===\n\n";
    fichier.close();
    afficher_population(population);

}
void remplacement(Population& population, const Population& enfants, int elite_count) {
    ofstream fichier("log.txt", ios::app);
    
    fichier << "\n=== DEBUT REMPLACEMENT ===\n";
    fichier << "Taille population actuelle: " << population.size() << "\n";
    fichier << "Nombre d'enfants: " << enfants.size() << "\n";
    fichier << "Nombre d'élites à conserver: " << elite_count << "\n";
    
       // Trie la population actuelle par fitness décroissante
    sort(population.begin(), population.end(), [](const Chromosome& a, const Chromosome& b) {
        return fonctionFitness(a) > fonctionFitness(b);
    });
        // Loggue les meilleurs individus actuels
    fichier << "\nMeilleurs individus actuels:\n";
    for (int i = 0; i < min(3, (int)population.size()); ++i) {
        fichier << i+1 << ". [ ";
        for (int g : population[i]) fichier << g << " ";
        fichier << "] Fitness: " << fonctionFitness(population[i]) << "\n";
    }
        // Trie les enfants par fitness décroissante
    Population enfants_tries = enfants;
    sort(enfants_tries.begin(), enfants_tries.end(), [](const Chromosome& a, const Chromosome& b) {
        return fonctionFitness(a) > fonctionFitness(b);
    });
        // Loggue les meilleurs enfants
    fichier << "\nMeilleurs enfants:\n";
    for (int i = 0; i < min(3, (int)enfants_tries.size()); ++i) {
        fichier << i+1 << ". [ ";
        for (int g : enfants_tries[i]) fichier << g << " ";
        fichier << "] Fitness: " << fonctionFitness(enfants_tries[i]) << "\n";
    }
    
        // Initialise la nouvelle population
    Population nouvelle_population;
    
       // Élitisme: conserve les 'elite_count' meilleurs individus
    for (int i = 0; i < elite_count && i < population.size(); ++i) {
        nouvelle_population.push_back(population[i]);
    }
    fichier << "\nÉlites conservées: " << nouvelle_population.size() << "\n";
        // Ajoute les meilleurs enfants jusqu'à remplir la population
    int index_enfant = 0;
    while (nouvelle_population.size() < population.size() &&
       static_cast<size_t>(index_enfant) < enfants_tries.size())
{
        nouvelle_population.push_back(enfants_tries[index_enfant]);
        index_enfant++;
    }
    fichier << "Enfants ajoutés: " << index_enfant << "\n";
        // Si nécessaire, complète avec de nouveaux individus aléatoires
    // Complément avec nouveaux individus si nécessaire
    if (nouvelle_population.size() < population.size()) {
        int manquants = population.size() - nouvelle_population.size();
        fichier << "Génération de " << manquants << " nouveaux individus\n";
        
        for (int i = 0; i < manquants; ++i) {
            Chromosome nouveau;
            for (int j = 0; j < N_CAPTEURS; ++j) {
                nouveau.push_back(rand() % N_EMPLACEMENTS);
            }
          //  reparer(nouveau);
            nouvelle_population.push_back(nouveau);
            
            fichier << "Nouvel individu " << i+1 << ": [ ";
            for (int g : nouveau) fichier << g << " ";
            fichier << "] Fitness: " << fonctionFitness(nouveau) << "\n";
        }
    }
        // Remplace l'ancienne population par la nouvelle
    population = nouvelle_population;
        // Calcule et loggue les statistiques de fitness
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

// Fonction pour réinitialiser les variables globales entre chaque instance
void reinitialiserVariables() {
    capteurs.clear();
    points_interet.clear();
    emplacements.clear();
    matrice_distance.clear();
    meilleure_solution.clear();
    meilleure_fitness = 0.0;
    TAILLE_POPULATION = 0;
    N_GENERATIONS = 0;
    MUTATION_RATE = 0.0;
    CROSSOVER_RATE = 0.0;
    TOURNAMENT_SIZE = 0;
    N_CAPTEURS = 0;
    N_EMPLACEMENTS = 0;
    chromosome_size = 0;
}

// Initialisation fichier CSV (sans écraser si déjà existant)
void initialiserFichierResultats(const string& nomFichier) {
    ifstream f_existant(nomFichier);
    if (!f_existant.good() || f_existant.peek() == ifstream::traits_type::eof()) {
        // Si le fichier n'existe pas ou est vide, on écrit l'entête
        ofstream f(nomFichier);
        f << "Instance,Meilleure_solution,POI_couverts,Nb_POI,Capteurs_utilises,Temps_ms\n";
        f.close();
    }
}

int main() {
    // Initialiser le fichier CSV
    initialiserFichierResultats("resultats.csv");
    // TEST: écriture forcée d'une ligne dans le CSV pour diagnostic
    {
        ofstream test_csv("resultats.csv", ios::app);
        if (test_csv.is_open()) {//Vérifie si le fichier est accessible
            //test_csv << "TEST,1,2,Oui,3,4,5,6\n";
            test_csv.close();//Ferme le fichier (écriture vide, juste pour tester)
            cout << "[DEBUG] Ligne de test écrite dans resultats.csv" << endl;
        } else {//Affiche une erreur si l'ouverture échoue
            cerr << "[ERREUR] Impossible d'ouvrir resultats.csv pour écriture de test." << endl;
        }
    }
    // Lire le fichier JSON  Ouvre le fichier JSON "data1.json" en lecture
    ifstream input("data1.json");
    if (!input.is_open()) {//Ouvre le fichier JSON "data1.json" en lecture
        cerr << "Erreur: Impossible d'ouvrir le fichier data1.json" << endl;//affiche une erreur et quitte le programme avec un code 1 (échec)
        return 1;
    }
    json j;//Crée un objet JSON pour stocker les données parsées
    try {
        input >> j;//Parse le contenu du fichier JSON dans l'objet `j`
        cout << "[DEBUG] Lecture JSON OK" << endl;//Confirme la réussite de la lecture
    } catch (const exception& e) {
        cerr << "Erreur lors de la lecture du fichier JSON: " << e.what() << endl;//Capture les erreurs de parsing et affiche le message
        return 1;
    }//Vérifie si la clé "genetic_params" existe dans le JSON
    // Récupérer les paramètres génétiques
    if (!j.contains("genetic_params")) {
        cerr << "Erreur: Pas de clé genetic_params dans le JSON." << endl;//Affiche une erreur si la clé est manquante
        return 1;
    }//Récupère les paramètres génétiques depuis l'objet JSON
    json genetic_params = j["genetic_params"];
    // Parcourir toutes les clés .bmp
    int instance = 0;
    for (auto it = j.begin(); it != j.end(); ++it) {
        std::string key = it.key();
        if (key.size() > 4 && key.substr(key.size()-4) == ".bmp") {
            cout << "\n=== Traitement de l'image " << key << " ===" << endl;
            reinitialiserVariables();//Réinitialise les variables globales pour une nouvelle instance
            auto start_time = chrono::high_resolution_clock::now();//Démarre un chronomètre pour mesurer le temps d'exécution
            try {
                lireInstanceDepuisObjet(it.value(), genetic_params);//Charge les données de l'instance (capteurs, POI, etc.)
                calculerMatriceDistance();
                calculerEtAfficherMatriceCouvertureParCapteur();
                Population population = initialiser_population();
                afficher_population(population);
               // Initialise les variables de contrôle de la convergence
                int generation = 0;
                int stagnation_count = 0;
                const int MAX_GENERATIONS = N_GENERATIONS;
                const int MAX_STAGNATION = 10;
                double meilleure_fitness_precedente = -1.0;
                Chromosome meilleure_solution_globale;
                double meilleure_fitness_globale = -1.0;
                while (generation < MAX_GENERATIONS && stagnation_count < MAX_STAGNATION) {//Boucle principale de l'algorithme génétique
                    log("\n--- Génération " + to_string(generation + 1) + " ---");
                    int solutions_valides = 0;
                    vector<double> fitness_values;
                    for (const auto& individu : population) {
                        //Calcule la fitness de l'individu
                        double fitness = fonctionFitness(individu);
                        fitness_values.push_back(fitness);
                       // Compte les solutions valides (fitness > 0)
                        if (fitness > 0) solutions_valides++;
                        if (fitness > meilleure_fitness_globale && estCouvertureComplete(individu)) {//Met à jour la meilleure solution globale
                            meilleure_fitness_globale = fitness;
                            meilleure_solution_globale = individu;
                        }
                    }//Affiche l'avancement en console
                    cout << "Génération " << generation + 1 
                         << " - Solutions valides: " << solutions_valides 
                         << "/" << population.size() << endl;
                    Population selected_population = selectionParTournoiAvecLog(population, TOURNAMENT_SIZE);
                    Population crossed_population = crossover(selected_population, CROSSOVER_RATE, chromosome_size);
                    mutation(crossed_population, MUTATION_RATE, N_EMPLACEMENTS);
                    for (auto& individu : crossed_population) {
                        if (!estCouvertureComplete(individu)) {
                            //reparer(individu);
                        }
                    }
                    remplacement(population, crossed_population, 2);
                    //Vérifie la stagnation (pas d'amélioration)
                    double max_fitness = *max_element(fitness_values.begin(), fitness_values.end());
                    if (max_fitness <= meilleure_fitness_precedente) {
                        stagnation_count++;
                        log("Stagnation #" + to_string(stagnation_count));
                    } else {
                        stagnation_count = 0;
                        meilleure_fitness_precedente = max_fitness;
                    }
                    generation++;//Passe à la génération suivante
                }//Arrête le chronomètre et calcule la durée
                auto end_time = chrono::high_resolution_clock::now();
                auto duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count();
                //Ouvre le fichier CSV en mode append
                ofstream result_csv("resultats.csv", ios::app);
                if (result_csv.is_open()) {
                    int nb_POI = points_interet.size();
                    int capteurs_utilises = 0;
                    int poi_couverts = 0;
                    string meilleure_config = "";
                    if (!meilleure_solution_globale.empty()) {//Si une solution valide a été trouvée Si une solution valide a été trouvée
                        capteurs_utilises = count_if(meilleure_solution_globale.begin(), meilleure_solution_globale.end(), [](int x) { return x != -1; });
                        poi_couverts = CountCoveredPOI(meilleure_solution_globale);//Compte les POI couverts
                        // Encoder la configuration sous forme de string
                        for (size_t i = 0; i < meilleure_solution_globale.size(); ++i) {//ormate la configuration pour le CSV
                            if (i > 0) meilleure_config += ";";
                            meilleure_config += to_string(meilleure_solution_globale[i]);
                        }
                    } else {
                        meilleure_config = "Aucune";
                    }//Écrit une ligne dans le CSV
                    result_csv << "[" << key << "],"
                              << "[" << meilleure_config << "],"
                              << "[" << poi_couverts << "],"
                              << "[" << nb_POI << "],"
                              << "[" << capteurs_utilises << "],"
                              << "[" << duration << "]\n";
                    result_csv.close();
                    // Affichage diagnostic Affiche des infos de débogage en console
                    cout << "[DIAG] Image: " << key << endl;
                    cout << "[DIAG] Nombre de points d'intérêt: " << nb_POI << endl;
                    cout << "[DIAG] Nombre de points couverts par la meilleure solution: " << poi_couverts << endl;
                    double rayon_max = 0.0;
                    for (const auto& c : capteurs) rayon_max = max(rayon_max, c.rayon);
                    cout << "[DIAG] Rayon maximal des capteurs: " << rayon_max << endl;
                    cout << "[DIAG] Configuration meilleure solution: ";
                    if (!meilleure_solution_globale.empty()) {
                        for (int v : meilleure_solution_globale) cout << v << " ";
                        cout << endl;
                    } else {
                        cout << "Aucune solution" << endl;
                    }
                }
            } catch (const exception& e) {// Capture les erreurs pendant le traitement d'une instance
                cerr << "[ERREUR] Exception dans le traitement de l'image " << key << ": " << e.what() << endl;
                ofstream result_csv("resultats.csv", ios::app);
                if (result_csv.is_open()) {
                    result_csv << "[" << key << "],"
                              << "[0],"
                              << "[0],"
                              << "[0],"
                              << "[0],"
                              << "[0]"
                              << "\n";
                    result_csv.close();
                }
                continue;//Passe à l'instance suivante
            }
            instance++;//Incrémente le compteur d'instances
        }
    }//Retourne 0 (succès) à la fin du programme.
    return 0;
}

