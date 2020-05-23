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
- `void ip_mat_free(ip_mat *a)`: Libera memoria se la matrice `a` è non nulla.
- `void ip_mat_init_random(ip_mat *t, float mean, float var)`: Genera una matrice di valori casuali, presi dalla matrice iniziale `t`, utilizzando la distribuzione normale di media `mean` e deviazione standard `std`.
- `ip_mat *ip_mat_sum(ip_mat *a, ip_mat *b)`: Esegue la somma due matrici se entrambe non sono nulle.
- `ip_mat *ip_mat_sub(ip_mat *a, ip_mat *b)`: Esegue la sottrazione due matrici se entrambe non sono nulle.
- `ip_mat *ip_mat_mean(ip_mat *a, ip_mat *b)`: Esegue la media tra due matrici, calcolando la media di ogni singolo elemento `v = (a+b)/2`, se entrambe le matrici non sono nulle.
- `ip_mat *ip_mat_mul_scalar(ip_mat *a, float c)`: Esegue la moltiplicazione scalare tra gli elementi di `a` e `c`, a patto che a non sia nulla.
- `ip_mat *ip_mat_add_scalar(ip_mat *a, float c)`: Esegue la somma scalare tra gli elementi di `a` e `c`, a patto che a non sia nulla.


### Parte 2

- `ip_mat *ip_mat_brighten(ip_mat *a, float bright)`: Aumenta o diminuisce la luminosità dell'immagine sommando ad ogni elemento di `a`, con `a` non nullo, il valore di `bright` e successivamente limitando i valori in `[0.0f,255.0f]`.
- `ip_mat *ip_mat_blend(ip_mat *a, ip_mat *b, float alpha)`: Fonde due immagini insieme (a patto che le loro matrici corrispondenti siano della stessa dimensione e non nulle), di una certa percentuale `alpha`, seguendo l'equazione `v = a*alpha + b*(1-alpha)` e riscalando i valori in modo tale che siano comprensi in `[0.0f, 255.0f]`.

### Parte 3

- `ip_mat *ip_mat_padding(ip_mat *a, int pad_h, int pad_w)`: Applica del padding all'immagine sui 4 lati seguendo l'equazione `dim = dim + 2*pad`, dove `dim` indica `a->h` o `a->w`.
- `void clamp(ip_mat *t, float low, float high)`: Limita i valori della matrice `t` non nulla tra `[low, high]`.
- `void rescale(ip_mat *t, float new_max)`: Scala i valori della matrice `t` non nulla prima in `[0.0f, 1.0f]`, poi in `[0, new_max]`.
- `float calculate_convolution(ip_mat *a, ip_mat *ker, int i, int j, int k)`: Calcola il valore di un elemento della matrice, con `a` e `ker` non nulli, seguendo la formula della convoluzione `val += val_a * val_ker`, dove `val_a` corrisponde al valore della matrice `a`, mentre `val_ker` al valore della matrice `ker` nella stessa posizione. `a` altro non è che una sottomatrice della matrice dell'immagine di dimensioni `ker->w * ker->h * ker->k`.
- `ip_mat *ip_mat_convolve(ip_mat *a, ip_mat *f)`: Esegue la convoluzione sulla matrice-immagine `a`, applicando il filtro `f`, entrambi non nulli.
- `ip_mat *create_average_filter(int w, int h, int k)`: Crea un filtro *average* di dimensioni `w*h*k`, seguendo la formula `val = 1/(w*h)`, il cui scopo è quello di creare un blur "fondendo" insieme i pixel adiacenti.
- `ip_mat *create_edge_filter()`: Crea un filtro `3x3x3` il cui scopo è quello di mostrare il contorno dell'immagine

### Extra
### Ausiliarie

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