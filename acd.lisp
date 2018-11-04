;;;; Decomposes arbitrary input polygon into convex sub-pieces using the ACD algorithm
(uiop:define-package #:cook/geojson/acd
    (:use #:cl
          #:alexandria
          #:iterate)
  (:export #:get-decomposition)
  (:documentation "Decomposes arbitrary polygons into convex sub-pieces"))
(in-package #:cook/geojson/acd)


;; appends point to the end of a list of points 
(defun add (point points-list)
  (append points-list (list point)))


;; return the index associated with the minimum element
(defun argmin (l)
  (let* ((min (apply #'min l)))
    (position min l :test #'=)))


;; return the index associated with the maximum element
(defun argmax (list)
  (let* ((max (apply #'max list)))
    (position max list :test #'=)))


;; accumulate all the x-coordinates
(defun get-x-coords (polygon)
  (mapcar #'first polygon))


;; accumulate all the y-coordinates
(defun get-y-coords (polygon)
  (mapcar #'second polygon))


;; return the euclidian distance between two points 
(defun get-distance (point1 point2)
  (sqrt (+ (expt (- (first  point2) (first  point1)) 2) 
           (expt (- (second point2) (second point1)) 2))))


;; return true if point is found within the hull list
(defun point-in-hull (point hull-list)
  (member point hull-list :test #'equal))


;; get the adjustment angle between previous point and next point
(defun get-adjustment-angle (heading previous-point next-point)
  (let* ((deg-in-circle 360)
         (left 270)
         (right 90)
         (delta-x (- (first  next-point) (first  previous-point)))
         (delta-y (- (second next-point) (second previous-point)))
         (relative-angle (* (/ 180 pi) (atan delta-y delta-x)))
         (adjustment-angle))

    ;; force relative angle to be measured as clockwise from East
    (if (< relative-angle 0) 
        (setf relative-angle (* -1 relative-angle)) 
        (setf relative-angle (- deg-in-circle relative-angle)))

    ;; force relative angle to be measured as clockwise from North
    (if (< relative-angle left) 
        (setf relative-angle (+ relative-angle right)) 
        (setf relative-angle (- relative-angle left)))

    ;; calculate adjustment angle
    (setf adjustment-angle (- relative-angle heading))
    (when (< adjustment-angle 0) 
      (setf adjustment-angle (+ adjustment-angle deg-in-circle))) 

    (values adjustment-angle)))


;; get the convex hull of the input polygon
(defun get-convex-hull (polygon)
  (let* ((left-point-index (argmin (get-x-coords polygon)))
         (left-point (nth left-point-index polygon))
         (hull-complete nil)
         (hull-list)
         (previous-point (nth left-point-index polygon))
         (heading 0)
         (current-adjustment)
         (current-distance)
         (current-point)
         (new-adjustment))
    
    ;; add the left-most point to the hull
    ;; TODO avoid constructing hull-list with `add`
    (setf hull-list (add left-point hull-list))

    ;; iterate until the convex hull is found
    (loop :while (not hull-complete)

          ;; initialize current values before new points are considered
          do  (setf current-adjustment 361)
              (setf current-distance 0)

              ;; iterate through the points to determine the next hull-point
              (loop :for next-point :in polygon

                    ;; consider all possible next points
                    do  (when (not (equal next-point previous-point))

                          ;; calculate the new adjustment angle for this point
                          (setf new-adjustment
                                (get-adjustment-angle heading previous-point next-point))

                          ;; maintain the best candidate for the next point in the hull list 
                          (when (or (< new-adjustment current-adjustment) 
                                    (and (= new-adjustment current-adjustment) 
                                         (< (get-distance previous-point next-point)
                                            current-distance)))
                            (setf current-adjustment  new-adjustment)
                            (setf current-point next-point)
                            (setf current-distance (get-distance previous-point next-point)))))

              ;; assess the hull complete condition
              (if (equal current-point left-point)
                  (setf hull-complete t)

                  ;; if not complete, add a new point to the hull
                  (progn
                    (setf hull-list (add current-point hull-list))
                    (setf heading (+ heading current-adjustment))
                    (setf previous-point current-point))))

    (values hull-list)))


;; get the concavity measurement of the input polygon
(defun get-concavity (polygon)  
  (let* ((hull (get-convex-hull polygon))
         (concavity-list)
         (left-bridge-index (length hull))
         (right-bridge-index 0)
         (max-concavity-index)
         (max-concavity)
         (notch))

    ;; TODO build concavity-list without using `add`
    ;; iterate through the points of the input polygon
    (loop :for point :in polygon
          do  (if (not (point-in-hull point hull))

                  ;; if point is not part of the hull
                  ;; TODO update how concavity is computed
                  (setf concavity-list (add (min (get-distance point (nth left-bridge-index hull))
                                                 (get-distance point (nth right-bridge-index hull)))
                                             concavity-list))

                  ;; if point is part of the hull
                  (progn
                    (setf concavity-list (add 0 concavity-list))

                    ;; adjust the bridge indices 
                    (setf right-bridge-index (+ right-bridge-index 1))
                    (if (= left-bridge-index (length hull))
                        (setf left-bridge-index 0)
                        (setf left-bridge-index (+ left-bridge-index 1))))))

    ;; extract the values to be returned
    (setf max-concavity-index (argmax concavity-list))
    (setf max-concavity       (nth max-concavity-index concavity-list))
    (setf notch               (nth max-concavity-index polygon))

    (values notch max-concavity concavity-list hull)))


;; cut the input polygon into two sub-pieces by adding a diagonal between notch and cut-point
(defun get-cut-polygon (polygon notch cut-point) 
  (let* ((index1 (position notch     polygon :test #'equal))
         (index2 (position cut-point polygon :test #'equal))
         (min-index (min index1 index2))
         (max-index (max index1 index2)))

         ;; extract polygon1 and polygon2 given the indices of the added diagonal
         (iter
           (for i :in (iota (length polygon)))
           (cond
             ((< i min-index)
              (collect (nth i polygon) :into polygon1))
             
             ((= i min-index)
              (collect (nth i polygon) :into polygon1)
              (collect (nth i polygon) :into polygon2))
             
             ((and (> i min-index) (< i max-index))
              (collect (nth i polygon) :into polygon2))

             ((= i max-index)
              (collect (nth i polygon) :into polygon1)
              (collect (nth i polygon) :into polygon2))

             ((> i max-index)
              (collect (nth i polygon) :into polygon1)))

           (finally
            (return-from get-cut-polygon (values (list polygon1 polygon2)))))))


;; get the interior angle at location given by middle point
(defun get-interior-angle (previous-point middle-point next-point)
  (let* ((deg-in-circle 360)
         (left 270)
         (reverse 180)
         (right 90)
         (interior-angle)
         (delta-x (- (first  middle-point) (first  previous-point)))
         (delta-y (- (second middle-point) (second previous-point)))
         (heading (* (/ 180 pi) (atan delta-y delta-x))))
    
    ;; force heading to be measured as clockwise from East
    (if (< heading 0)
        (setf heading (* heading -1))
        (setf heading (- deg-in-circle heading)))

    ;; force heading to be measured as clockwise from North
    (if (< heading left)
        (setf heading (+ heading right))
        (setf heading (- heading left)))

    ;; calculate the interior angle at the notch point (used to assess resolved condition)
    (setf interior-angle (- reverse (get-adjustment-angle heading middle-point next-point)))
    (values interior-angle)))


;; determine if notch r has been successfully resolved in every polygon in polygons
(defun valid-resolve (polygons notch)
  (let* ((previous-index)
         (next-index)
         (interior-angle)
         (num-sides-triangle 3))

    ;; return false if one of the polygons has less than three sides
    (unless (< (min (length (first polygons)) (length (second polygons))) num-sides-triangle)
      
       ;; iterate over the vertices of the polygon
      (loop :for polygon :in polygons
            :for notch-index = (position notch polygon :test #'equal)
            :do 
               
               ;; get the previous index before the notch point
               (setf previous-index (- notch-index 1))
               (when (= previous-index -1)
                 (setf previous-index (- (length polygon) 1)))

               ;; get the next point after the notch point
               (setf next-index (+ notch-index 1))
               (when (= next-index (length polygon))
                 (setf next-index 0))

               ;; get the interior angle at the notch point
               (setf interior-angle (get-interior-angle (nth previous-index polygon) 
                                                        (nth notch-index    polygon)
                                                        (nth next-index     polygon)))

            never (< interior-angle 0)))))


;; add an additional diagonal line such that notch is resolved and the score function is maximized
(defun get-resolved-polygons (input-polygon notch concavity-list)
  (let* ((Sc 0.1)
         (Sd 1.0)
         (score)
         (polygons)
         (best-score 0)
         (best-polygons))

    ;; iterate through the points of the input polygon to find the best diagonal to add
    (loop :for i :in (iota (length input-polygon))
          do  (when (not (equal notch (nth i input-polygon)))

                ;; get the score of this point
                (setf score (/ (+ 1 (* Sc (nth i concavity-list)))
                               (+ Sd (get-distance notch (nth i input-polygon)))))

                ;; if the new score is better, get the resulting polygons of adding this diagonal
                (when (> score best-score)
                  (setf polygons (get-cut-polygon input-polygon notch (nth i input-polygon)))

                  ;; when the resolve is valid, maintain the current best diagonal
                  (when (valid-resolve polygons notch)
                    (setf best-polygons (get-cut-polygon input-polygon notch (nth i input-polygon)))
                    (setf best-score score)))))

    (values best-polygons)))


;; take the approximate convex decomposition of polygon p
(defun acd (input-polygon tolerance reference-length)
  (multiple-value-bind (notch concavity concavity-list convex-hull)
      (get-concavity input-polygon)
    
    (let* ((resolved-polygons)
           (new-polygons)
           (output-polygons '()))

      ;; return input polygon if concavity tolerance is already satisfied
      (when (<= concavity (* tolerance reference-length))
        (return-from acd (list convex-hull)))

      ;; resolve the input polygon 
      (setf resolved-polygons (get-resolved-polygons input-polygon notch concavity-list))

      ;; recursively call the ACD algorithm
      (loop :for resolved-polygon :in resolved-polygons
            :do (setf new-polygons (acd resolved-polygon tolerance reference-length))
               (loop :for new-polygon :in new-polygons
                     :do  (push new-polygon output-polygons)))

      (values output-polygons))))


;; decomposes the input polygon into smaller convex sub-pieces
;; TODO pre-process input to work with ACD function
(defun get-decomposition (input-polygon tau length)
  (list input-polygon tau length)
  ;; Notes on current implementation
  ;; + acd operates on 'open' polygons
  ;; + acd operates in CW manner
  ;; + acd outputs 'open' polygons

  ;; Open TODOs
  ;; TODO -> eliminated use of `append` and `add` function
  ;; TODO -> update how concavity is computed
  )


;; run test on ACD algorithm 
(defun test ()
  (let* ((polygons)
         (tolerance)
         (notch)
         (cut-point)
         (valid)
         (reference-length)
         (test-count 0)
         (test-pass  0)
         (polygon-toy (list
                       (list  1.0  1.0)
                       (list  2.0  4.0)
                       (list  1.0  6.0)
                       (list  2.0  8.0)
                       (list  5.0  7.0)
                       (list  6.0 10.0)
                       (list 10.0 10.0)
                       (list 12.0 12.0)
                       (list 13.0 16.0)
                       (list 13.0 18.0)
                       (list 16.0 18.0)
                       (list 16.0 16.0)
                       (list 15.0 10.0)
                       (list 16.0  7.0)
                       (list 16.0  4.0)
                       (list 14.0  2.0)
                       (list 12.0  3.0)
                       (list 13.0  5.0)
                       (list 12.0  5.0)
                       (list  8.0  3.0)
                       (list 12.0  1.0)
                       (list  8.0  2.0)
                       (list  2.0  1.0)))
         
         (polygon-3d (list
                      (list  1.0  1.0  2.0)
                      (list  2.0  4.0  5.5)
                      (list  1.0  6.0  6.6)
                      (list  2.0  8.0 10.0)
                      (list  5.0  7.0 12.0)
                      (list  6.0 10.0  0.0)
                      (list 10.0 10.0  0.0)
                      (list 12.0 12.0  8.0)
                      (list 13.0 16.0 10.0)
                      (list 13.0 18.0 64.0)
                      (list 16.0 18.0 64.9)
                      (list 16.0 16.0  9.9)
                      (list 15.0 10.0  0.8)
                      (list 16.0  7.0  0.0)
                      (list 16.0  4.0  0.0)
                      (list 14.0  2.0  0.0)
                      (list 12.0  3.0  4.0)
                      (list 13.0  5.0  3.3)
                      (list 12.0  5.0 10.0)
                      (list  8.0  3.0 11.1)
                      (list 12.0  1.0 12.0)
                      (list  8.0  2.0  5.0)
                      (list  2.0  1.0  0.0))))

    ;; Test 1: decomposition of simple polygons, tolerance = 0
    (setf tolerance 0.0)
    (setf reference-length 16.0)
    (setf polygons (acd polygon-toy tolerance reference-length))
    (assert (= 9 (length polygons)))
    (loop :for p :in polygons
          :do (assert (>= (length p) 3)))
    (format t ">> Test 1: Passed")
    (setf test-count (+ 1 test-count))
    (setf test-pass  (+ 1 test-pass))

    
    ;; Test 2: decomposition of simple polygons, tolerance = 1
    (setf tolerance 1.0)
    (setf reference-length 16.0)
    (setf polygons (acd polygon-toy tolerance reference-length))
    (assert (= 1 (length polygons)))
    (loop :for p :in polygons
          :do (assert (>= (length p) 3)))
    (format t "~%>> Test 2: Passed")
    (setf test-count (+ 1 test-count))
    (setf test-pass  (+ 1 test-pass))
    

    ;; Test 3: tests `cut-polygon` and `valid-resolve` functions
    (setf notch      (list 10.0 10.0))
    (setf cut-point  (list 15.0 10.0))
    (setf polygons   (get-cut-polygon polygon-toy notch cut-point))
    (setf valid      (valid-resolve polygons notch))
    (assert valid)
    (setf test-count (+ 1 test-count))
    (setf test-pass  (+ 1 test-pass))
    (format t "~%>> Test 3: Passed")

    
    ;; Test 4: decomposition on simple 3D polygon, tolerance = 0
    (setf tolerance 0.0)
    (setf reference-length 16.0)
    (setf polygons (acd polygon-3d tolerance reference-length))
    (assert (= 9 (length polygons)))
    (loop :for p :in polygons
          :do (assert (>= (length p) 3)))
    (format t "~%>> Test 4: Passed")
    (setf test-count (+ 1 test-count))
    (setf test-pass  (+ 1 test-pass))

    
    ;; Test 5: decomposition on simple 3D polygon, tolerance = 0.25
    (setf tolerance 0.25)
    (setf reference-length 16.0)
    (setf polygons (acd polygon-3d tolerance reference-length))
    (assert (= 6 (length polygons)))
    (loop :for p :in polygons
          :do (assert (>= (length p) 3)))
    (format t "~%>> Test 5: Passed")
    (setf test-count (+ 1 test-count))
    (setf test-pass  (+ 1 test-pass))

    (format t "~%----------~%>> ~S of ~S tests are passing" test-pass test-count)))







