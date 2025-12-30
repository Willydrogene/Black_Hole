#include <stdio.h>
#include <stdlib.h> 
#include <math.h>
#include <omp.h>    // Permet le calcul parallèle



//////////////////////////////////////////////// Fonction qui permet d'afficher une barre de progression ////////////////////////////////////////////////
void progress_bar(double progress)
{
    const int barWidth = 50;                // Longueur de la barre (nombre de caractères)
    
    if (progress < 0.0) progress = 0.0;     // On s'assure que progress reste dans [0, 1]
    if (progress > 1.0) progress = 1.0;

    int pos = (int)(progress * barWidth);   // Position actuelle dans la barre (index de la "tête" de la progression)

    printf("\r[");                          // '\r' ramène le curseur au début de la même ligne pour réécrire par-dessus
 

    for (int i = 0; i < barWidth; i++) {    // On parcourt chaque "case" de la barre
        if (i < pos)        printf("=");    // Partie déjà complétée : on affiche des '='
        else if (i == pos)  printf(">");    // Position actuelle : on met un '>' pour symboliser l'avancement
        else                printf(" ");    // Partie restante : on laisse des espaces vides
    }

    printf("] %6.2f%%", progress * 100.0);  // Affiche le pourcentage avec 2 décimales, sur 6 caractères de large
    fflush(stdout);                         // On force l'affichage immédiat (évite que ça reste dans le buffer)
}



//////////////////////////////////////////////////////////////// Equation des géodésiques ////////////////////////////////////////////////////////////////
double du_point_dtheta(double u) {  // u'' = 3/2 u^2 - u 
	return ((3.0/2.0) * u * u) - u; // Notre équation de trajectoire
}



/////////////////////// Fonction qui permet de calculer theta_d en fonction de alpha et i, c'est l'équation qui est dans le sujet ///////////////////////
double theta_d(double alpha,double i){
    if( 1.0 - cos(alpha)*cos(alpha)*cos(i)*cos(i) <= 1e-12) {     // On vérifie simplement si le dénominateur est égal à 0
        return acos( (-sin(alpha)*cos(i)) / sqrt(1e-12) );
    }
    return acos( (-sin(alpha)*cos(i)) / sqrt(1.0-cos(alpha)*cos(alpha)*cos(i)*cos(i)) );  // Donné dans le sujet
}



////////////////////////////////////////////////////////////////// Runge Kutta ordre 4 //////////////////////////////////////////////////////////////////
double RK4(double *u, double *u_point, double dtheta, int Ntheta, double rs) {
    double r=rs/ *u;                                        // On calcule ce qui nous intéresse pour le futur plot

	for (int i=0; i<Ntheta; i++) {                          // Boucle sur tous les thétas
		double k1 = (*u_point) * dtheta;                            // Premier terme u de RK4
		double k1_point = du_point_dtheta(*u) * dtheta;             // Premier terme u' de RK4

		double k2 = (*u_point + k1_point/2.0) * dtheta;             // Deuxième terme u de RK4
		double k2_point = du_point_dtheta(*u + k1/2.0) * dtheta;    // Deuxième terme u' de RK4

		double k3 = (*u_point + k2_point/2.0) * dtheta;             // Troisième terme u de RK4
		double k3_point = du_point_dtheta(*u + k2/2.0) * dtheta;    // Troisième terme u' de RK4

		double k4 = (*u_point + k3_point) * dtheta;                 // Quatrième terme u de RK4
		double k4_point = du_point_dtheta(*u + k3) * dtheta;        // Quatrième terme u' de RK4
		
		*u = *u + (k1 + 2*k2 + 2*k3 + k4)/6.0;                                        // u final de RK4
		*u_point = *u_point + (k1_point + 2*k2_point + 2*k3_point + k4_point)/6.0;    // u' final de RK4

        r=rs/ *u;                         // On calcule r à la fin de notre tour de boucle                                            
        if (r<rs || r>1e8) break;       // On regarde si on est tombé dans le trou noir ou si on est très très loin pour couper le calcul 
	}
    return r;   // On retourne r pour l'exploiter en dehors de la fonctoin
}



////////////////////////////////////////////////////// Equation de luminosité du disque d'accrétion //////////////////////////////////////////////////////
double LumDisque(double r, double rs) {
    const double c = 299792458; // Vitesse de la lumière
    const double m_point = 2.0; // Débit d'accrétion de masse

    double x = r / (3.0 * rs);  // Coordonnée radiale adimensionnée 
    if (x <= 1.0) return 0.0;   // Vérification du domaine de validité (pour éviter le NaN/Inf)
    
    double F0 = (3.0 / 4.0) * (m_point * c * c) / (4.0 * M_PI * rs * rs);       // Calcul de F0 : (3 * m_point * c^2) / (4 * pi * rs^2)
    double ratio = ( 1.0 / ((x*x*x) * (1.0 - (1.0 / (2.0 * x)) )) )     
                   * (1.0 - (1.0 / sqrt(x)) 
                   + ( sqrt(1.0 / (8.0 * x)) * log( ((sqrt(2.0 * x) + 1.0) / (sqrt(2.0 * x) - 1.0)) * ((sqrt(2.0) - 1.0) / (sqrt(2.0) + 1.0)) ) ) );  // Calcul du ratio
    
    double F_emis = F0 * ratio;     // Calcul de F_emis
    return F_emis;                  // On retourne la aleur de F_emis
}



///////////////////////////////////////////////////////////////////// Equation de 1+z /////////////////////////////////////////////////////////////////////
double UnplusZ(double uc, double i, double alpha, double r, double rs) {
    double x = r / (3.0 * rs);          // Coordonnée radiale adimensionnée 
    if (x <= 1.0) return 1.0;           // Le flux est nul, mais pour éviter la division par zéro dans F_obs, on retourne 1.

    double T1 = 1.0 / sqrt(1.0 - (1.0 / (2.0 * x)));        // Terme 1 (Gravitational Redshift/Dilatation du temps)
    double T2_magnitude = 1.0 / (sqrt(8.0*x*x*x));          // Terme 2 de magnitude
    double T2_angular = (1.0 / uc) * cos(i) * cos(alpha);   // Terme 2 angulaire, inclut la vitesse et les angles
    
    return T1 * (1.0 + T2_magnitude * T2_angular);          // On retourne la valeur de 1 + z = ...
}



////////////////////////////////////////////////////////////// Calculs des images secondaires //////////////////////////////////////////////////////////////
void images_n_disque(int nb_images, FILE* file, double bc, double b, double alpha, double i, double theta_start, double dtheta, double rs, double R_min, double R_max) {
    double uc = bc/b;                                           // Variable sans dimension donnée dans le sujet
    double u = 1e-3;  						    				// r0 = rs/u0; lumière qui arrive de loin, r0 -> inf, donc u0 << 1 ici on prend par ex 1e-3
    double u_point = sqrt( (4.0/27.0) * uc*uc - u*u + u*u*u);   // Première valeur de u', à partir de l'équation adimensionnée
    double theta_end = theta_d(alpha,i);                        //condition d'arrêt de l'intégration pour le photon considéré lorsqu'il touche le disque
    int Ntheta = ((theta_end)-theta_start)/dtheta;;	            // Nombre d'élément théta à calculer
    int a_touche_le_disque = 0;                                 // Savoir si on a touché le disque, on l'initialise à 0 par défaut (comme un false)
    
    for (int n=0; n<nb_images; n++){                                                                        // On boucle sur le nombre d'images qu'on veut (secondaire, tertiaire etc...)
        if(a_touche_le_disque == 0){                                                        			    // Si le disque n'a pas été touché (false)

            if (n>0) Ntheta = ((theta_end + (n * M_PI))-theta_end)/dtheta;                                  // Nombre d'élément à calculer (+n*pi ici pour voir si on va toucher le disque plus loin)
            
            double r = RK4(&u, &u_point, dtheta, Ntheta, rs);                          			            // On récupère la dernière valeur de r après avoir résolu l'equa diff 

            if(r>R_min && r<R_max){                                                         			    // Condition pour savoir si on est sur le disque d'accrétion
                double F_emis = LumDisque(r,rs);                                                   			// Valeur du Flux émis par le disque
                double redshift = UnplusZ(uc,i,alpha,r,rs);                                        			// Valeur de 1+z
                double redshift2 = UnplusZ(uc,i,-alpha+M_PI,r,rs);                                        	// Valeur de 1+z pour la partie symétrique
                double F_obs = F_emis/(redshift * redshift * redshift * redshift);                 			// Valeur de F observé, F_obs = F_emis/(1+z)^4
                double F_obs2 = F_emis/(redshift2 * redshift2 * redshift2 * redshift2);                 	// Valeur de F observé pour la partie symétrique
                
                #pragma omp critical(file_write)                                                            // Section critique: un seul thread à la fois peut écrire dans le fichier
                {
                    fprintf(file,"%.5lf \t %.5lf \t %.20e \t %.20e \t %.20e \n", b, alpha, r, F_obs, redshift);         // On écrit toutes les valeurs dans le fichier
                    fprintf(file,"%.5lf \t %.5lf \t %.20e \t %.20e \t %.20e \n", b, -alpha+M_PI, r, F_obs2, redshift2); // On écrit toutes les valeurs symétriques dans le fichier
                }
                a_touche_le_disque = 1;                                                         			// On a touché le disque alors on met la valeur 1 (true)
            }
        }
    }
}



////////////////////////////////////////////////////////////////////// Fonction Exploration des Conditions Initiales //////////////////////////////////////////////////////////////////////
void explorationCI(double i, double b_end, double b_start, double db, double alpha_start, double daplha, double theta_start, double dtheta, double rs, double R_min, double R_max, int nb_images) {

    FILE* file = fopen("image_detecteur.txt","w");                                              // On ouvre un fichier .txt pour écrire dedans
    fprintf(file, "\t %-10s %-20s %-30s %-30s %-30s\n", "b", "alpha", "r", "F_obs", "UnpZ");    // Entête de notre fichier .txt

    const double bc = 3.0 * sqrt(3.0) * rs /2.0;        // constante sans dimension donnée dans le sujet
    
    int n_steps = (int)((b_end - b_start) / db);        // Pour la barre de progression, calcule le nombre total 
    
    int progress_counter = 0;                           // Compteur pour la barre
    
    #pragma omp parallel                                // Début de la région parallèle : les threads OpenMP sont lancés ici
    {
        #pragma omp for schedule(dynamic)                       // La boucle sur 'step' est répartie entre les threads, avec un ordonnancement dynamique
        for (int step = 0; step < n_steps; ++step) {                    // step est privé à chaque thread, géré par OpenMP
            double b = b_start + step * db;                             // On reconstruit b à partir de step
            
            #pragma omp atomic                                          // Incrément atomique: évite les conditions de course sur 'progress_counter'
            progress_counter++;
            
            #pragma omp critical(progress_print)                        // Affichage de la barre pour une seule thread à la fois
            {
                progress_bar(progress_counter / (double)n_steps);                 // Barre de progression qui se met à jour
            }
            
            for(double alpha = alpha_start ; alpha<3.0*M_PI/2.0 ; alpha+=daplha){                               // Boucle sur alpha (demi image), (séquentielle à l'intérieur de chaque thread)                     
                images_n_disque(nb_images, file, bc, b, alpha, i, theta_start, dtheta, rs, R_min, R_max);       // Fonction qui calcule les nb_images du trou noir
            }
        }
    }
    printf("\nTerminé !\n");    // Tous les threads ont fini leurs tâches
    fclose(file);               // On n'oublie pas de fermer le fichier une fois terminé
}



////////////////////////////////////////////////////////////////////// MAIN //////////////////////////////////////////////////////////////////////
int main(){
    double start = omp_get_wtime();         // Début chronométrage (temps réel)

    const double rs = 1;                    // Rayon de Schwarzschild

    const double i = M_PI/50;               // Angle i, inclinaison du disque d'accrétion

    const double b_start = 0.01;            // Première valeur de b du détecteur
    const double b_end = 25*rs;             // Dernière valeur de b (on arrête b à cette valeur), c'est la taille du détecteur
    const double db = 1e-2;                 // Le pas de b, ce qui va déterminer le nombre d'élément b, + il est petit, + les calculs sont longs
    
    const double alpha_start = M_PI/2.0;    // Angle alpha du détecteur, première valeur de alpha (on ne fait qu'une demie image)
    const double daplha = 1e-2;             // Le pas de alpha, idem que db
    
    const double theta_start = 0.001;       // Première valeur de l'angle theta
    const double dtheta = 1e-3;             // Le pas de theta, idem que db et dalpha
    
    const double R_min = 3*rs;              // Rayon du disque d'accrétion minimum qu'on fixe à 3rs pour avoir trou noir -> disque autour du trou noir
    const double R_max = 20*rs;             // Rayon maximal du disque d'accrétion

    const int nb_images = 4;                // Nombre d'images secondaires du trou noir

    explorationCI(i, b_end, b_start, db, alpha_start, daplha, theta_start, dtheta, rs, R_min, R_max, nb_images);  // Appel de la fonction qui fait les calculs
    
    
    double end = omp_get_wtime();                               // Fin chronométrage (temps réel)
    double elapsed = end - start;                               // Calcul du temps mis (en secondes)
    printf("Temps total : %.3f secondes\n", elapsed);           // Printf du temps total
    #pragma omp parallel
    {
        #pragma omp single                                      // Une seule thread affiche l’information sur le nombre de threads utilisés
        {
            printf("Threads utilisés : %d\n", omp_get_num_threads());
        }
    }

    return 0;                               // On retourne 0 à la fin du main
}
