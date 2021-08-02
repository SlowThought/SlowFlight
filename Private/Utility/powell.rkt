#lang racket/base
(require "./1d-search.rkt" ; for powell-1, fv->fr
         math/number-theory
         racket/future
         racket/math
         racket/vector)
(provide powell-1); powell)
(module+ test
  (require plot "../test-private.rkt")
  (printf"powell.rkt: tests running.~n"))

; Powell too much algorithm. Just throw CPU cycles at the problem
;#| powell-n finds minimum of real function of vector of reals of arbitrary dimension.
;   Initially, search is along normal unit vectors. Search results are used to
;   modify the search vectors to approximate a gradient descent.
;|#
;
;(define (powell-n f x0 h ϵ max)
;  ; f(x) is function to be minimized - x is a vector
;  ; x0 is starting position, a vector
;  ; h is initial step size along search vector
;  ; epsilon is step size limit, if we move less than epsilon in a given direction, we stop
;  ; max is the max number of iterations of the algorithm
;  (define n (vector-length x0))
;  (define n+1 (add1 n))
;  ; Search vectors
;  (define e (build-basis-vectors n))
;  (define ss e)
;  ; Search history
;  (define xs (make-vector n+1 (make-vector n 0)))
;  (vector-set! xs 0 x0)
;  (define fs (make-vector n+1))
;  (vector-set! fs 0 (f x0))
;  ; Big loop
;  (for [(m (in-range max))]
;    (cond [(= (remainder m n))
;           (set! ss e)])
;    ; Perform search
;    (for*[(i (in-range max))
;          (j (in-range n))]
;      (define k (add1 j))
;      (let-values([(x* f*)
;                   (search-along-direction (vector-ref ss j)
;                                           (vector-ref xs j)
;                                           f h ϵ max)])
;        (vector-set! xs k x*)
;        (vector-set! fs k f*)))
;    ; Identify least fruitful search direction
;    (let [(k 0) ;k for kill this vector
;          (dxk (vector-diff^2 (vector-ref xs 0)(vector-ref xs 1)))]
;      (for [(i (in-range 1 n))]
;        (let [(dxi (vector-diff^2 (vector-ref xs i)(vector-ref xs (add1 i))))]
;          (cond [(< dxi dxk)
;                 (set! k i)
;                 (set! dxk dxi)])))
;      ; Replace least fruitful with first
;      (vector-set! ss k (vector-ref ss 0)))
;    ; Replace first with best (from x0 to xn-2, we exclude last search direction because best should be normal to it)
;    (vector-set! ss 0 (let*[(dx (vector-map - (vector-ref xs (- n 2)) (vector-ref xs 0)))
;                            (ndx (sqrt (vector-diff^2 (vector-ref xs (- n 2)) (vector-ref xs 0))))]
;                        (vector-map (λ (e) (/ e ndx)) dx)))
;    ; Repeat search with new search vectors; replace x0, f0
;    (vector-set! xs 0 (vector-ref xs n))
;    (vector-set! fs 0 (vector-ref fs n)))
;  (printf "powell-n likely isn't done. Returns recent search results for now.~n")
;  (values xs fs))

;; powell helper functions

; build-basis-vectors generates n normal vectors of dimension n
(define (build-basis-vectors n)
  (build-vector n (λ (i)
                    (define v (make-vector n 0.))
                    (vector-set! v i 1.)
                    v)))
  
; measure magnitude of difference between vectors as sum of squares of differences of elements
(define (vector-diff^2 v1 v2)
  (for/fold [(return 0.)]
            [(x1 (in-vector v1))
             (x2 (in-vector v2))]
    (+ return (*(- x1 x2)(- x1 x2)))))
  
(module+ test
  (check-equal? (build-basis-vectors 2) #(#(1.0 0.0)#(0.0 1.0)))
  (check-equal? (vector-diff^2 #(1 2 3) #(2 3 2)) 3.))

; test infrastructure -- generate random nd quadratic equations
(module+ test
  (require math/distributions math/matrix)
  ; Generate a random positive definite matrix
  (define (random-pd-matrix n)
    (define nd (normal-dist 0. 0.5))
    (define A
      (build-matrix    
       n n
       (λ (i j) (cond [(= i j)
                       (exp (sample nd))] ; Diagonal must be positive
                      [(> i j) 0]         ; Building triangular matrix to ensure diagonal determines determinant
                      [else (sample nd)]))))
    ;(printf "A : ~s~n" A) ; uncomment to debug. confidence is high.
    (matrix+ A (matrix-transpose A))); Suspect, but have not proven, that this last step, making final matrix symmetrical, is unneccesary.
  ; Generate a random quadratic function that has a minimum f* at location x*
  (define (gen-random-f f* x*)
    (define n (vector-length x*))
    (define A (random-pd-matrix n))
    (define X* (->col-matrix x*))
    (λ (x)
      (define X (->col-matrix x))
      (define dX (matrix- X X*))
      (+ f* (matrix-ref (matrix* (matrix-transpose dX) (matrix* A dX))
                       0 0))))
  ; Tests of gen-random-f, 2D examples consistently show expected result
    (let*[(f (gen-random-f 1. #(1. 1.)))
          (g (λ (x y) (f (vector x y))))]
      (display (plot3d (surface3d g)
                       #:x-min 0
                       #:x-max 2
                       #:y-min 0
                       #:y-max 2))
      (printf "Above you should see a surface with a clear minimum near (1,1)~n")
      (for-each (λ (x)
                   (printf "~nX: ~a f(X): ~a " x (f x)))
                (list #(1. 1.) #(1. 2.) #(2. 1.) #(1. 0.) #(0. 1.))))
  (printf "If the above values of x and f(x) make sense, then ~nbelieve the random test cases are working.~n"))

; search-along-direction applies powell-1 to a function of a vector
(define (search-along-direction s x f h ϵ max)
  (let*[(g (fv->fr f x s))
        (dt (powell-1 g (- h) 0 h ϵ max))
        (x* (vector-map + x (vector-map (λ(e)(* dt e)) s)))
        (f* (f x*))]
    (values x* f*)))

(module+ test
  ; Test search-along-direction
  (let[(f (λ (v) (let [(x (vector-ref v 1))] (* (sub1 x)(sub1 x)))))];only second element is significant, minimum at v=#(0 1 0)(
    (define-values (x* f*)(search-along-direction #(0 1 0) #(0 0 0) f 0.5 0.01 5))
    (check-= f* 0 .01)
    (check-= (vector-ref x* 0) 0 .01)
    (check-= (vector-ref x* 1) 1 .01)
    (check-= (vector-ref x* 2) 0 .01)))

#| Experimental algorithms. "king-s" (for sequential) is meant to be equivalent to one cycle of powell-n. "king-p" (for parallel)
   is meant to use the same pattern but exploit parallelism for a speedup roughly proportional to the number of processors available
   (assuming dimension of problem is larger than number of processors).

   The King algorithm, for both sequential and parallel, is

   - Search along each basis vector
   - Search a "best" vector that is normal to last/best search
   - Search the last/best direction again

   Unlike powell, the search vectors are not updated with conjugate vectors. Rather, we take the two final "good" search directions
   as good enough, and then start over. Less book keeping should compensate for lost efficiency (if it exists) in original Powell.

|#

(define sad-max 10) ; Max iterations of search-along-direction

(define (king-s f x0 ss h ϵ n)
  ; f is function to be minimized
  ; x0 is starting point of search
  ; ss are search directions
  ; h is an initial step size
  ; ϵ is a tolerance, used as a stopping criterion
  ; n is the dimension of x0 and ss
  ;; Search each direction in ss sequentially
  (define x x0)
  (let [(xs (for/list [(s ss)]
            (let-values ([(x* f*)(search-along-direction s x f h ϵ sad-max)])
              (set! x x*)
              x*)))]
    ; Search best direction (approx gradient excluding last search)
    (define best-s (vector-map - (cadr xs) x0))
    (let-values ([(x* f*)(search-along-direction best-s x f h ϵ sad-max)])
      ; Search last direction (only available normal to best)
      (let-values ([(x* f*)(search-along-direction (vector-ref ss (sub1 n)) x* f h ϵ sad-max)])
        ; Return best input, minimized output
        (values x* f*)))))

(module+ test
  (printf "Sequential search test - true min function 3.14, true min location #(1 2 3 4)~n")
  (let [(f (gen-random-f 3.14 #(1. 2. 3. 4.)))
        (x0 #(0. 0. 0. 0.))
        (e (build-basis-vectors 4))
        (ϵ 0.01)
        (h 0.1)
        (n 4)]
    (for [(i (in-range 5))]
      (let-values ([(x* f*)(king-s f x0 e h ϵ n)])
        (printf "X=~a, f(X)=~a~n" x* f*)
        (set! x0 x*)))
    (printf "Best X:~a, best f(X):~a~n" x0 (f x0))))
; 21-07-23 above tested good.

(define (king-p f x0 ss h ϵ n)
  ; f is function to be minimized
  ; x0 is starting point of search
  ; ss are search directions
  ; h is an initial step size
  ; ϵ is a tolerance, used as a stopping criterion
  ; n is the dimension of x0 and ss
  
  (define f0 (f x0))
  ; Search each direction in ss in parallel
  (define futures (vector-map (λ (s) (future (λ () (search-along-direction s x0 f h ϵ sad-max)))) ss))
  ; Collect results
  (define xs (make-vector n))
  (define fs (make-vector n))
  (for [(i (in-range n))]
    (let-values ([(x f)(touch (vector-ref futures i))])
      (vector-set! xs i x)
      (vector-set! fs i f)))
  ; Find best direction
  (define i* 0)
  (define f* (vector-ref fs 0))
  (for [(i (in-range 1 n))]
    (cond [(< (vector-ref fs i) f*)
           (set! f* (vector-ref fs i))
           (set! i* i)]))
  ; Start at best result of parallel search, exclude that direction
  (define s* (build-vector n (λ (i) (if (= i i*) 0.
                                        (- (vector-ref (vector-ref xs i) i)
                                           (vector-ref x0 i))))))
  (set!-values (x0 f*)(search-along-direction s* (vector-ref xs i*) f h ϵ sad-max))
  ; One more search along best direction
  (search-along-direction (vector-ref ss i*) x0 f h ϵ sad-max))

(module+ test
  (printf "Parallel search test - early days.~n")
  (let [(f (gen-random-f 3.14 #(1. 2. 3. 4.)))
        (x0 #(0. 0. 0. 0.))
        (e (build-basis-vectors 4))
        (ϵ 0.01)
        (h 0.1)]
    (for [(i (in-range 5))]
      (let-values ([(x* f*)(king-p f x0 e h ϵ 4)])
        (printf "X=~a, f(X)=~a~n" x* f*)
        (set! x0 x*)))
    (printf "Best X:~a, best f(X):~a~n" x0 (f x0)))
  (printf "~nParallel seems to work almost as well as sequential, but is it faster?~n")
  (define x0 #(0. 0. 0. 0.))
  (define e (build-basis-vectors 4))
  (define ϵ 0.01)
  (define h 0.1)
  (for[(king-x (list king-s king-p))
       (algorithm (list "Sequential" "Parallel"))]
    (printf "~s timing results:~s~n"
            algorithm
            (time (for [(i (in-range 1000))]
                    (let [(f(gen-random-f 3.14 #(1. 2. 3. 4.)))]
                      (king-x f x0 e ϵ h 4)))))))
; 21-08-01 Parallel is no faster - are we achieving parallelism at all?
           
  
  
