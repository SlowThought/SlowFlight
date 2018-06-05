#lang racket/base
(provide integrate)
;#| f is a linear function between x1 & x2 (f = mx + b).
;   Find definite integral of f^n in the interval [x1, x2]. |#
;
;(define(integrate f1 f2 x1 x2 n)
;  (define dx (- x2 x1))
;  (if (= dx 0.) 0.
;      (let [(m (/ (- f2 f1) dx))]
;        (if (= m 0) 
;            (* (expt f1 n) dx)
;            (let [(n+1 (add1 n))]
;              (/ (- (expt f2 n+1)(expt f1 n+1))
;                 m n+1))))))
; 
;(module+ test
;  (require "test-private.rkt")
;  (printf "quad.rkt running tests.~n")
;  ; delta x = 0
;  (check-= (integrate 1. 2. 3. 3. 7) 0. epsilon)
;  ; f = 1, F = x
;  (check-= (integrate 1. 1. 2. 3. 1) 1. epsilon)
;  (check-= (integrate 2. 2. 0. 1. 2) 4. epsilon)
;  (check-= (integrate 3. 3. 6. 7. 3) 27. epsilon)
;  ; f = x + b, F = 1/2 x^2 + bx
;  (check-= (integrate 0. 1. 0. 1. 1) .5 epsilon)
;  (check-= (integrate 0. 2. 0. 2. 1) 2. epsilon)
;  (check-= (integrate 0. 2. 2. 4. 1) 2. epsilon)
;  (check-= (integrate 0. 1. 0. 1. 2) (/ 3.) epsilon)
;  (check-= (integrate 0. 1. 7. 8. 2) (/ 3.) epsilon))

#| NOTE I AM BREAKING EVERY USER OF INTEGRATE -- but I've broken everything else.
   "integrate" finds definite integral of f and its derivatives over interval h
   by evaluating the Taylor series (TS) (whats the plural of series?)
   f(x+h) = f(x) + f'(x)h + f''(x)h^2/2 + f'''(x)h^3/6 + ...
   f'(x+h) = f'(x)+ f''(x)h + ...
   f''(x+h) = f''(x) + f'''(x)h + ...
   ...

   A succinct representation: f(x+h)=Σf^(i)(x) * h^i/i!

   f^(i)(x) represents the i-th derivative of f at x.

|#

(define (integrate h . fs)
  (if (= (length fs) 1) fs ; The last (highest) derivative is assumed constant over the interval [x, x+h]
                           ; We ignore (null? fs)
      ; We assume (> (length fs) 1)
      (cons ; We return a list of this derivative's value at the next step and all the higher derivatives at the next step
       ; Evaluate TS of this derivative; sum f^(i)(x) * h^i/i!
       (let loop [(sum 0.)
                  (i+1 1) ; i == 0
                  (h^i 1.); h^0 == 1
                  (i! 1)
                  (fi (car fs))
                  (fr (cdr fs))]
         (if (null? fr) (+ sum (/(* fi h^i)i!))
             (loop (+ sum (/(* fi h^i)i!))
                   (add1 i+1)
                   (* h^i h)
                   (* i! i+1)
                   (car fr)
                   (cdr fr))))
       ; Eval TS of higher derivatives
       (apply integrate (cons h (cdr fs))))))

                     
(module+ test
  (require plot racket/math "../test-private.rkt")
  (printf "quad.rkt running tests.~n")
  (printf "y = x~n")
  (plot (list (function (λ (x) x))
              (points (let [(h 0.1)
                            (f 0.0)
                            (fp 1.)]
                        (for/list[(i (in-range 10))]
                          (let [(fs (integrate h f fp))]
                          (set! f (car fs))
                          (vector (* h (add1 i))f))))))
        #:x-min 0.
        #:x-max 1.)
  ; f = x^2; f' = 2x; f" = 2; x0 = 0. f(h)= .01; f'(h)=.02; f" = 2
  (printf "y = x^2~n")
  (plot (list (function (λ (x)(* x x)))
              (points (let [(h   .1)
                            (f   0.)
                            (fp  0.)
                            (fpp 2.)]
                        (for/list[(i (in-range 10))]
                          (let [(fs (integrate h f fp fpp))]
                            (set! f (car fs))
                            (set! fp(cadr fs))
                            (vector (* h (add1 i))f))))))
        #:x-min 0.
        #:x-max 1.)
  ; f = sin x; f' = cos x; f" = - sin x;
  (printf "y = sin x~n")
  (plot (cons (function sin)
              (let* [(n 20)
                     (h  (/(* 2 pi) n))]
                (list (points (let [(f   0.)
                                    (fp  1.)]
                        (for/list[(i (in-range n))]
                          (let [(fs (integrate h f fp))]
                            (set! f (car fs))
                            (set! fp(cos (* h (add1 i))))
                            (vector (* h (add1 i))f)))))
                      (points (let [(f   0.)
                                    (fp  1.)
                                    (fpp 0.)]
                                (for/list[(i (in-range n))]
                                  (let [(fs (integrate h f fp fpp))]
                                    (set! f (car fs))
                                    (set! fp(cadr fs))
                                    (set! fpp (- f))
                                    (vector (* h (add1 i))f))))
                              #:sym 'triangle))))                        
        #:x-min 0.
        #:x-max (* 2 pi)))