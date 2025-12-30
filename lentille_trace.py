import numpy as np
import matplotlib.pyplot as plt


# ========================= PARAMÈTRES =========================
BIN_FILE = "image_detecteur.bin"    # Fichier binaire généré par le C (float32, 3 valeurs par record)
SOURCE_IMAGE = "fond.png"           # Image de fond (doit être dans le même dossier)
OUT_PNG = "image_lentille.png"      # Sortie enregistrée

b_max = 70/np.sqrt(2)               # La taille du détecteur dans le C est à 70*rs format circulaire, mais en divisant par racine(2) ici ça permet l'image en carré

det_res = 5000                      # Résolution de l'image du détecteur (plus élevé = moins d'artefacts, mais nécessite plus de rayons)


xbg_min, xbg_max = -100.0, 100.0      # Fenêtre physique du plan de fond z = -D, qui correspond à l'étendue couverte par l'image SOURCE_IMAGE
ybg_min, ybg_max = -100.0, 100.0

# Composition des pixels
# "average" : moyenne si collisions (souvent plus joli)
# "overwrite": dernier rayon écrase (plus "puriste" mais peut strier)
COMPOSITING = "average"             # Si det_res élevé, on voit quasi pas la différence entre les deux, s'il est très bas, la différence se voit

# Lecture streaming du binaire
# Un record = 3 float32 = 12 octets.
# chunk_records = nb de records lus par itération.
chunk_records = 3_000_000           # ~36 Mo de data brute par chunk (safe dans WSL)



# ========================= CHARGEMENT IMAGE SOURCE =========================
src = plt.imread(SOURCE_IMAGE)                      # Lecture de l'image représentant le plan de fond (z = -D).


if src.dtype == np.uint8:                           # On convertit systématiquement l'image en float32 normalisé dans [0, 1]
    src_f = src.astype(np.float32) / 255.0
else:
    src_f = src.astype(np.float32)
    if src_f.max() > 1.0:                           # Par sécurité : certaines images float peuvent être codées sur [0, 255]
        src_f /= 255.0

# --- grayscale -> RGB
if src_f.ndim == 2:                                 # Si l'image est 2D (H × W), on la convertit en RGB (H × W × 3)
    src_f = np.stack([src_f, src_f, src_f], axis=2)

H, W, C = src_f.shape                               # H : hauteur de l'image, W : largeur de l'image, C : nombre de canaux (3 pour RGB, 4 pour RGBA)


# ========================= BUFFERS DETECTEUR =========================
if COMPOSITING == "average":                                                # Mode "average" : On accumule la somme des contributions lumineuses reçues par chaque pixel du détecteur, ainsi que le nombre de rayons qui y arrivent.
    acc = np.zeros((det_res, det_res, C), dtype=np.float64)                 # acc[j, i, :] contient la somme des couleurs reçues par le pixel (i, j)
    cnt = np.zeros((det_res, det_res), dtype=np.int32)                      # cnt[j, i] contient le nombre de rayons ayant contribué au pixel (i, j)
elif COMPOSITING == "overwrite":                                            
    out = np.zeros((det_res, det_res, C), dtype=np.float32)                 # Mode "overwrite" : Chaque rayon écrase directement la valeur du pixel correspondant sur le détecteur (pas de moyenne en cas de collisions).
else:
    raise ValueError("COMPOSITING doit être 'average' ou 'overwrite'.")     # Sécurité : empêche l'exécution si une méthode de composition inconnue est fournie


# ========================= OUTILS DE MAPPING =========================
def coord_to_index(x, xmin, xmax, n):
    t = (x - xmin) / (xmax - xmin)              # Normalisation de la coordonnée x dans l'intervalle [0, 1]
    return (t * (n - 1)).astype(np.int32)       # Conversion de la coordonnée normalisée en indice de pixel. Le facteur (n - 1) assure que xmin -> 0 et xmax -> n - 1


# ========================= LECTURE STREAMING + RENDU =========================
floats_per_chunk = chunk_records * 3    # Chaque rayon (record) est stocké sous la forme de 3 float32 : (b, alpha, y_local). On lit donc 3 * chunk_records flottants par itération.
total_records = 0                       # Compteur du nombre total de rayons traités (utile pour suivre l'avancement).


with open(BIN_FILE, "rb") as f:         # Ouverture du fichier binaire en lecture ("rb" = read binary). Le streaming permet de traiter des fichiers très volumineux sans les charger en RAM.
    
    while True:
        buf = np.fromfile(f, dtype=np.float32, count=floats_per_chunk)  # Lecture d'un bloc de float32 depuis le fichier. np.fromfile lit directement le binaire (très rapide, sans parsing texte).
        if buf.size == 0:
            break                                                       # Fin de fichier : aucune donnée lue.
        if buf.size % 3 != 0:
            buf = buf[: (buf.size // 3) * 3]                            # Par sécurité : si la dernière lecture ne contient pas un nombre entier de records (fichier tronqué ou fin incomplète), on coupe l'excédent.

        data = buf.reshape(-1, 3)               # Reformatage en tableau (N, 3) : N rayons, 3 valeurs par rayon.
        b = data[:, 0]
        alpha = data[:, 1]
        y_local = data[:, 2]                    # coordonnée transverse au point d'intersection avec le plan de fond (z = -D) (mesurée dans le plan 2D du rayon, puis réorientée par alpha)

        # --- coordonnées détecteur
        x_det = b * np.cos(alpha)               # Le détecteur est échantillonné en coordonnées polaires (b, alpha)
        y_det = b * np.sin(alpha)               # On reconstruit la position cartésienne (x_det, y_det) associée au pixel du détecteur.

        # --- coordonnées sur le plan du fond
        x_bg = y_local * np.cos(alpha)          # Le code C intègre la géodésique dans un plan 2D (plan du rayon).
        y_bg = y_local * np.sin(alpha)          # On "tourne" ce plan dans le repère global (X,Y) grâce à l'angle azimutal alpha.

        # --- pixels détecteur
        i_det = coord_to_index(x_det, -b_max, b_max, det_res)   # Mapping linéaire : x_det,y_det ∈ [-b_max, b_max] -> indices ∈ [0, det_res-1].
        j_det = coord_to_index(y_det, -b_max, b_max, det_res)   # (i_det = colonne, j_det = ligne)

        # --- pixels source (dans l'image de fond)
        u_src = coord_to_index(x_bg, xbg_min, xbg_max, W)       # Mapping linéaire : (x_bg,y_bg) dans la fenêtre physique du fond -> pixels de SOURCE_IMAGE.
        v_src = coord_to_index(y_bg, ybg_min, ybg_max, H)       # u_src = colonne dans l'image source, v_src = ligne dans l'image source.

        # --- masque bornes
        mask = (                                # Ce masque élimine :
            (i_det >= 0) & (i_det < det_res) &  # - les rayons dont l'impact détecteur est en dehors de la zone reconstruite
            (j_det >= 0) & (j_det < det_res) &
            (u_src >= 0) & (u_src < W) &        # - les rayons qui pointent hors de la portion du fond représentée par SOURCE_IMAGE
            (v_src >= 0) & (v_src < H)
        )

        # --- Application du masque sur tous les tableaux d'indices pour rester cohérent
        i_det = i_det[mask]
        j_det = j_det[mask]
        u_src = u_src[mask]
        v_src = v_src[mask]

        # --- récupération des couleurs sur l'image source
        colors = src_f[v_src, u_src, :]                         # Pour chaque rayon restant : on lit la couleur du pixel touché sur le plan de fond.

        # --- accumulation des contributions sur le détecteur
        if COMPOSITING == "average":
            np.add.at(acc, (j_det, i_det, slice(None)), colors) # np.add.at gère correctement les collisions : plusieurs rayons peuvent tomber
            np.add.at(cnt, (j_det, i_det), 1)                   # sur le même pixel du détecteur. On somme les couleurs et on compte les contributions.
        else:
            out[j_det, i_det, :] = colors                       # Mode overwrite : le dernier rayon écrit la couleur du pixel (pas de moyenne).

        # --- suivi d'avancement
        total_records += data.shape[0]                          # data.shape[0] = nombre de records lus dans ce chunk (avant masque)
        if total_records % (30_000_000) < data.shape[0]:        # Petit affichage périodique pour suivre la progression sur gros fichiers
            print(f"Traité ~{total_records:,} records")


# ========================= FINALISATION =========================
if COMPOSITING == "average":
    out = np.zeros((det_res, det_res, C), dtype=np.float32)                 # Image finale du détecteur (float32, valeurs attendues dans [0,1])
    nonzero = cnt > 0                                                       # Masque des pixels effectivement touchés par au moins un rayon
    out[nonzero] = (acc[nonzero] / cnt[nonzero, None]).astype(np.float32)   # Moyenne des contributions : (somme des couleurs) / (nombre de contributions) 
                                                                            # cnt[nonzero, None] ajoute un axe pour diffuser correctement la division sur les canaux couleur.


out = np.clip(out, 0.0, 1.0)                                                # (évite tout dépassement dû à des erreurs d'arrondi ou à des images mal normalisées)


# ========================= AFFICHAGE + SAUVEGARDE =========================
plt.figure(figsize=(8, 8))                          # Crée une figure carrée pour visualiser l'image reconstruite au détecteur
plt.imshow(out)                                     # Affiche l'image
plt.axis("off")                                     # Masque les axes pour obtenir une visualisation "pure image"
plt.title("Image lentillée (binaire streaming)")    # Titre de la figure (optionnel)
#plt.show()                                          # Affiche la figure à l'écran

# --- Sauvegarde PNG 
out8 = (out * 255.0 + 0.5).astype(np.uint8)         # Ici on convertit explicitement en uint8 pour une sauvegarde robuste.
plt.imsave(OUT_PNG, out8)                           # Sauvegarde dans le dossier

# --- Infos en plus
print(f"Enregistré : {OUT_PNG}")
print(f"Records traités : {total_records:,}")
print(f"det_res={det_res}, b_max={b_max}, bg_window=[{xbg_min},{xbg_max}]x[{ybg_min},{ybg_max}], COMPOSITING={COMPOSITING}")
