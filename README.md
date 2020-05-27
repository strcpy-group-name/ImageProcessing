# ImageProcessing

## Gruppo
Gruppo: strcpy(group, "name");

## Descrizione

### Parte 1
- `float get_val(ip_mat *a, unsigned int i, unsigned int j, unsigned int k)`:
Restituisce il valore della matrice `a`, contenuto in posizione (i,j,k) se la matrice non è nulla e gli indici sono compresi in `[0,dim]`, dove `dim` è rispettivamente `a->h,a->w,a->k`.
- `void set_val(ip_mat *a, unsigned int i, unsigned int j, unsigned int k, float v)`:
Pone il valore della matrice `a` in posizione (i,j,k) a `v` se la matrice non è nulla e gli indici sono compresi in `[0,dim]`, dove `dim` è rispettivamente `a->h,a->w,a->k`.
- `void compute_indexes(unsigned int i, unsigned int* ih, unsigned int* iw, unsigned int* ik, unsigned int w, unsigned int k)`: Calcola, dato un indice monodimensionale `i`, gli indici tridimensionali della matrice e il salva in `ih,iw,ik`.
- `ip_mat *ip_mat_create(unsigned int h, unsigned int w, unsigned int k, float v)`: Dato il numero di righe `h`, il numero di colonne `w` e il numero di canali `k` crea una matrice inizializzando tutti i valori della matrice a `v` e assegnando righe, colonne e canali ai parametri `h`, `w`, `k`. La matrice 3D è completamente linearizzata.
- `void ip_mat_free(ip_mat *a)`: Libera memoria se la matrice `a` è non nulla.
- `void ip_mat_init_random(ip_mat *t, float mean, float var)`: Genera una matrice di valori casuali, presi dalla matrice iniziale `t`, utilizzando la distribuzione normale di media `mean` e deviazione standard `std`.
- `ip_mat *ip_mat_subset(ip_mat *t, unsigned int row_start, unsigned int row_end, unsigned int col_start, unsigned int col_end)`:Data una matrice `t`restituisce la sottomatrice di `t` a partire dalla riga `row_start` e colonna `col_start` incluse fino alla riga `row_end` e colonna `col_end` escluse. 
- `ip_mat *ip_mat_sum(ip_mat *a, ip_mat *b)`: Esegue la somma due matrici se entrambe non sono nulle.
- `ip_mat *ip_mat_sub(ip_mat *a, ip_mat *b)`: Esegue la sottrazione due matrici se entrambe non sono nulle.
- `ip_mat *ip_mat_mean(ip_mat *a, ip_mat *b)`: Esegue la media tra due matrici, calcolando la media di ogni singolo elemento `v = (a+b)/2`, se entrambe le matrici non sono nulle.
- `ip_mat *ip_mat_mul_scalar(ip_mat *a, float c)`: Esegue la moltiplicazione scalare tra gli elementi di `a` e `c`, a patto che a non sia nulla.
- `ip_mat *ip_mat_add_scalar(ip_mat *a, float c)`: Esegue la somma scalare tra gli elementi di `a` e `c`, a patto che a non sia nulla.


### Parte 2

- `ip_mat *ip_mat_brighten(ip_mat *a, float bright)`: Aumenta o diminuisce la luminosità dell'immagine sommando ad ogni elemento di `a`, con `a` non nullo, il valore di `bright` e successivamente limitando i valori in `[0.0f,255.0f]`.
- `ip_mat *ip_mat_to_gray_scale(ip_mat *in)`: Data una matrice `in` non nulla, converte l'immagine in scala di grigi. In ogni pixel viene calcolata la media dei 3 canali e il valore ottenuto viene replicato su tutti e 3 i canali.
- `ip_mat *ip_mat_blend(ip_mat *a, ip_mat *b, float alpha)`: Fonde due immagini insieme (a patto che le loro matrici corrispondenti siano della stessa dimensione e non nulle), di una certa percentuale `alpha`, seguendo l'equazione `v = a*alpha + b*(1-alpha)` e riscalando i valori in modo tale che siano comprensi in `[0.0f, 255.0f]`.

### Parte 3

- `ip_mat *ip_mat_padding(ip_mat *a, int pad_h, int pad_w)`: Applica del padding all'immagine sui 4 lati seguendo l'equazione `dim = dim + 2*pad`, dove `dim` indica `a->h` o `a->w`.
- `void clamp(ip_mat *t, float low, float high)`: Limita i valori della matrice `t` non nulla tra `[low, high]`.
- `void rescale(ip_mat *t, float new_max)`: Scala i valori della matrice `t` non nulla prima in `[0.0f, 1.0f]`, poi in `[0, new_max]`.
- `float calculate_convolution(ip_mat *a, ip_mat *ker, int i, int j, int k)`: Calcola il valore di un elemento della matrice, con `a` e `ker` non nulli, seguendo la formula della convoluzione `val += val_a * val_ker`, dove `val_a` corrisponde al valore della matrice `a`, mentre `val_ker` al valore della matrice `ker` nella stessa posizione. `a` altro non è che una sottomatrice della matrice dell'immagine di dimensioni `ker->w * ker->h * ker->k`.
- `ip_mat *ip_mat_convolve(ip_mat *a, ip_mat *f)`: Esegue la convoluzione sulla matrice-immagine `a`, applicando il filtro `f`, entrambi non nulli.
- `ip_mat *create_gaussian_filter(int w, int h, int k, float sigma)`: Date le dimensioni del kernel e il valore `sigma`, restituisce il filtro gaussiano da sottomettere alla convolve. Si calcola la cella centrale del kernel avente coordinate `cx` e `cy` e per ciascuna posizione viene calcolato il valore da inserire secondo la formula: 1/(2πσ^2) * e^((x^2+y^2)/2σ^2), dove x e y sono la distanza dal centro del kernel e `σ` la costante passata come parametro di funzione che determina l'effetto della vicinanza dei pixel al centro del kernel.
- `ip_mat *create_average_filter(int w, int h, int k)`: Crea un filtro *average* di dimensioni `w*h*k`, seguendo la formula `val = 1/(w*h)`, il cui scopo è quello di creare un blur "fondendo" insieme i pixel adiacenti.
- `ip_mat *create_edge_filter()`: Crea un filtro `3x3x3` il cui scopo è quello di mostrare il contorno dell'immagine
- `ip_mat *create_emboss_filter()`: Crea un filtro `3x3x3` che dà all'immagine un effetto di profondità/£D.

### Extra
- `ip_mat *ip_mat_to_gray_scale_lum_corr(ip_mat *in)`: La scala di grigi non è più effettuata facendo la media dei 3 canali ma facendone una media pesata tenendo conto della diversa percezione dei colori da parte dell'occhio umano. Ai 3 canali viene dato un peso diverso. Questo approccio migliora la qualità dell'immagine e lo si può verificare facilmente soprattutto con immagini scure, che normalmente subirebbero un ulteriore scurimento che con questo metodo non avviene.
- `ip_mat *ip_mat_to_gray_scale_gamma_corr(ip_mat *in)`: le immagini normalmente sono gamma-compresse in uno spazio RGB non lineare, per effettuare la conversione è quindi necessario fare una gamma-estensione, calcolare la media pesata dei canali e poi rieseguire la gamma-compressione. Con questo metodo si ottiene un'immagine in scala di grigi molto fedele all'originale.

### Ausiliarie
- `void compute_indexes(unsigned int i, unsigned int* ih, unsigned int* iw, unsigned int* ik, unsigned int w, unsigned int k)` è una funzione ausiliaria per il calcolo degli indici necessaria nei cicli for, dove abbiamo usato un unico indice e serve una conversione rapida per l'uso di alcune funzioni elementari come set e get.
- `float calculate_convolution(ip_mat *a, ip_mat *ker, int i, int j, int k)` è una funzione ausiliaria di `ip_mat_convolve`.
- `float gamma_correction_exp_to_linear(float v)` è una funzione ausiliaria di `ip_mat_to_gray_scale_gamma_corr` ed esegue un'espansine gamma per convertire il singolo canale in uno spazio lineare.

## Funzioni extra
Ogni eventuale aggiunta extra aggiungetela in fondo a `ip_lib.c` sotto `/*** FUNZIONI EXTRA ***/` e in `ip_lib.h` fate lo stesso.

## Funzioni ausiliarie
In caso di funzioni ausiliarie aggiungere `/**** FUNZIONE AUSILIARIA ***/` sopra la descrizione della funzione in `ip_lib.c`.

## Usage

Per compilare tutto usare
```sh
make all
```
oppure 
```sh
make
```

Per compilare ed eseguire i test
```
make test
./test
```

Per compilare ed eseguire il programma
```sh
make main
./main <params>
```

Sono inclusi due script di test `runner_gauss.sh` e `runner_noise.sh`
Per eseguirli usare
```sh
make main
./runner_<gauss|noise>.sh
```

Per verificare eventuali memory leak
```sh
make debug
valgrind -v --leak-check=full ./debug
```