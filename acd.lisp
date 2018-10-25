;; Implementation of the Approximate Convex Decomposition (ACD) algorithm

;; utilize geojson.io for visualizing geojson files
;; eventual repository to be placed in: https://git.mers.csail.mit.edu/mars-toolkit/path-planning-models

;; get the range of indices from 0 to max, used for indexing lists
(defun range (max &key (min 0) (step 1))
  (loop for n from min below max by step
    collect n))


;; appends point (a list) to the end of list l 
(defun add (point l)
  (append l (list point)))


;; return the index associated with the minimum element
(defun argmin (l)
  (let* ((min (apply #'min l))
         (min-i))
  (loop for i in (range (length l))
    when (and (= min (nth i l))
              (not min-i))
      do (setf min-i i)
    finally (return  min-i))))


;; return the index associated with the minimum element
(defun argmax (l)
  (let* ((max (apply #'max l))
         (max-i))
  (loop for i in (range (length l))
    when (and (= max (nth i l))
              (not max-i))
      do (setf max-i i)
    finally (return  max-i))))


;; round a number to a specified decimal place
(defun round-to (number precision &optional (fun #'round))
  (let ((div (expt 10 precision)))
    (/ (funcall fun (* number div)) div)))


;; accumulate all the x-coordinates
(defun get-x-coords (p)
  (loop for xy in p
    collect (first xy)))


;; accumulate all the y-coordinates
(defun get-y-coords (p)
  (loop for xy in p
    collect (second xy)))


;; return the euclidian distance between two points 
(defun get-dist (point1 point2)
  (sqrt (+ (expt (- (first point2)  (first point1)) 2) 
           (expt (- (second point2) (second point1)) 2))))


;; return true iff point is found within the list hull
(defun point-in-hull (point hull)
  (loop for hull-point in hull
    do (when (and (= (first hull-point)  (first point))
                  (= (second hull-point) (second point)))
          (return t))
    finally (return nil)))


;; get the adjustment angle between previous point and next point
(defun get-adj-angle (heading previous-point next-point)
  (let* ((deg-in-circle 360)
         (left-turn 270)
         (right-turn 90)
         (sigfigs 5)
         (delta-x (- (first next-point)  (first previous-point)))
         (delta-y (- (second next-point) (second previous-point)))
         (rel-angle (round-to (* (/ 180 pi) (atan delta-y delta-x)) sigfigs))
         (angle-adj))

  ;; force rel-angle to be measured as clockwise from East
  (if (< rel-angle 0) 
    (setf rel-angle (* -1 rel-angle)) 
    (setf rel-angle (- deg-in-circle rel-angle)))

  ;; force rel-angle to be measured as clockwise from North
  (if (< rel-angle left-turn) 
    (setf rel-angle (+ rel-angle right-turn)) 
    (setf rel-angle (- rel-angle left-turn)))

  ;; calculate adjustment angle
  (setf angle-adj (- rel-angle heading))
  (when (< angle-adj 0) 
    (setf angle-adj (+ angle-adj deg-in-circle))) 

  (values angle-adj)))


;; get the convex hull of polygon p
(defun convex-hull (p)
  (let* ((hull-list)
         (hull-incomplete t)
         (left-point-index (argmin (get-x-coords p)))
         (left-point (list (nth left-point-index (get-x-coords p)) 
                           (nth left-point-index (get-y-coords p))))
         (previous-point (list (nth left-point-index (get-x-coords p)) 
                               (nth left-point-index (get-y-coords p))))
         (heading 0)
         (current-adj)
         (current-dist)
         (current-point)
         (new-adj))

  ;; add the left-most point to the hull
  (setf hull-list (add left-point hull-list))

  ;; iterate until the convex hull is found
  (loop while hull-incomplete

    ;; initialize current values before new points are considered
    do  (setf current-adj 360)
        (setf current-adj 360)
        (setf current-dist 0)

        ;; iterate through the points to determine the next hull-point
        (loop for next-point in p

          ;; if next-point != previous-point
          do  (when (not (and (= (first  next-point) (first previous-point)) 
                              (= (second next-point) (second previous-point))))

                  ;; calculate the new adjustment angle for this point
                  (setf new-adj (get-adj-angle heading previous-point next-point))

                  ;; maintain the best canidate for the next hull point 
                  (when (or (< new-adj current-adj) 
                            (and (= new-adj current-adj) 
                                 (< (get-dist previous-point next-point) current-dist)))
                      (setf current-adj new-adj)
                      (setf current-point (list (first next-point) (second next-point)))
                      (setf current-dist (get-dist previous-point next-point)))))

      ;; assess the hull complete condition
      (if (and (= (first left-point)  (first current-point)) 
               (= (second left-point) (second current-point)))
        (setf hull-incomplete nil)

        ;; add a new point to the hull
        (progn
          (setf hull-list (add current-point hull-list))
          (setf heading (+ heading current-adj))
          (setf previous-point (list (first current-point) (second current-point))))))

  (values hull-list)))


;; get the concavity measurement of polygon p
(defun concave (p)
  (let* ((hull (convex-hull p))
         (c-list)
         (left-bridge-index (length hull))
         (right-bridge-index 0)
         (max-c-index)
         (max-c)
         (notch))

  ;; iterate through the points of p
  (loop for point in p
    do  (if (not (point-in-hull point hull))

          ;; if point is not part of the hull
          (setf c-list (add (min (get-dist point (nth left-bridge-index hull))
                                 (get-dist point (nth right-bridge-index hull)))
                            c-list))

          ;; if point is part of the hull
          (progn
            (setf c-list (add 0 c-list))
            (setf right-bridge-index (+ right-bridge-index 1))
            (if (= left-bridge-index (length hull))
              (setf left-bridge-index 0)
              (setf left-bridge-index (+ left-bridge-index 1))))))

  ;; extract the values to be returned
  (setf max-c-index (argmax c-list))
  (setf max-c (nth max-c-index c-list))
  (setf notch (list (first  (nth max-c-index p))
                    (second (nth max-c-index p))))

  (values (list notch max-c c-list))))


;; return the new polygons created when diagonal between notch and point is added to p 
(defun get-resolved-polygons (p notch point)
  (let* ((index1)
         (index2)
         (polygon1)
         (polygon2))
  ;; extract the indices of notch and point
  (loop for i in (range (length p))
    do  (when (and (= (first notch)  (first (nth i p)))
                   (= (second notch) (second (nth i p))))
          (setf index1 i))
        (when (and (= (first point)  (first (nth i p)))
                   (= (second point) (second (nth i p))))
          (setf index2 i)))

  ;; extract polygon1 and polygon2 given the indices above
  (loop for i in (range (length p))
    do  (if (< i (min index1 index2))
          (setf polygon1 (add (list (first  (nth i p)) (second (nth i p))) polygon1))

          (if (= i (min index1 index2))
            (progn 
              (setf polygon1 (add (list (first  (nth i p)) (second (nth i p))) polygon1))
              (setf polygon2 (add (list (first  (nth i p)) (second (nth i p))) polygon2)))

            (if (and (> i (min index1 index2)) (< i (max index1 index2)))
              (setf polygon2 (add (list (first  (nth i p)) (second (nth i p))) polygon2))

              (if (= i (max index1 index2))
                (progn
                  (setf polygon1 (add (list (first  (nth i p)) (second (nth i p))) polygon1))
                  (setf polygon2 (add (list (first  (nth i p)) (second (nth i p))) polygon2)))

                (when (> i (max index1 index2))
                  (setf polygon1 (add (list (first  (nth i p)) (second (nth i p))) polygon1))))))))

  (values (list polygon1 polygon2))))


;; get the interior angle found at the notch point
(defun get-interior-angle (prev-point notch-point next-point)
  (let* ((deg-in-circle 360)
         (left-turn 270)
         (invert-angle 180)
         (right-turn 90)
         (interior-angle)
         (sigfigs 5)
         (delta-x (- (first notch-point)  (first prev-point)))
         (delta-y (- (second notch-point) (second prev-point)))
         (heading (round-to (* (/ 180 pi) (atan delta-y delta-x)) sigfigs)))

  ;; force heading to be measured as clockwise from East
  (if (< heading 0)
    (setf heading (* heading -1))
    (setf heading (- deg-in-circle heading)))

  ;; force heading to be measured as clockwise from North
  (if (< heading left-turn)
    (setf heading (+ heading right-turn))
    (setf heading (- heading left-turn)))

  ;; calculate the interior angle at the notch point (used to assess resolved condition)
  (setf interior-angle (- invert-angle (get-adj-angle heading notch-point next-point)))
  (values interior-angle)))


;; determine if notch r has been successfully resolved in every polygon in polygons
(defun valid-resolve (polygons notch)
  (let* ((valid t)
         (notch-index)
         (prev-index)
         (next-index)
         (interior-angle)
         (num-sides-triangle 3))

  ;; return false if one of the polygons has less than three sides
  (when (<  (min (length (first polygons)) (length (second polygons))) num-sides-triangle)
    (return-from valid-resolve nil))

  ;; iterate over polygons 
  (loop for polygon in polygons

    ;; iterate over the vertices of the polygon
    do  (loop for i in (range (length polygon))

          ;; extract the notch index for this polygon
          do  (when (and (= (first (nth i polygon))  (first notch))
                         (= (second (nth i polygon)) (second notch)))
                (setf notch-index i)))

        ;; get the previous index before the notch point
        (setf prev-index (- notch-index 1))
        (when (= prev-index -1)
          (setf prev-index (- (length polygon) 1)))

        ;; get the next point after the notch point
        (setf next-index (+ notch-index 1))
        (when (= next-index (length polygon))
          (setf next-index 0))

        ;; get the interior angle at the notch point to determine if it has been resolved properly
        (setf interior-angle (get-interior-angle (nth prev-index polygon) 
                                                 (nth notch-index polygon)
                                                 (nth next-index polygon)))

        ;; no longer valid of notch is not resolved in either polygon
        (when (< interior-angle 0)
          (setf valid nil)))

  (values valid)))


;; add an additional diagonal line such that notch r is resolved and the score function is maximized
;;  + returns the two new polygons achieved by splitting p with the optimal diagonal
(defun resolve (p notch c-list)
  (let* ((Sc 0.1)
         (Sd 1.0)
         (score)
         (polygons)
         (best-score 0)
         (best-polygons))

  ;; iterate through the points of p to find the best diagonal to add
  (loop for i in (range (length p))
    do  (when (not (and (= (first notch)  (first  (nth i p)))
                        (= (second notch) (second (nth i p)))))

          ;; get the score of this point
          (setf score (/ (+ 1 (* Sc (nth i c-list)))
                         (+ Sd (get-dist notch (nth i p)))))

          ;; if the new score is better, get the resulting polygons of adding this diagonal
          (when (> score best-score)
            (setf polygons (get-resolved-polygons p notch (nth i p)))

            ;; when the resolve is valid, maintain the current best diagonal
            (when (valid-resolve polygons notch)
              (setf best-polygons (get-resolved-polygons p notch (nth i p)))
              (setf best-score score)))))

  (values best-polygons)))


;; take the approximate convex decomposition of polygon p
(defun acd (p tau)
  (let* ((concave-return (concave p))
         (notch  (nth 0 concave-return))
         (c      (nth 1 concave-return))
         (c-list (nth 2 concave-return))
         (polygons)
         (polygon-list)
         (new-polygons))

  ;; return p if concavity tolerance is already satisfied
  (when (< c tau)
    (return-from acd (list p)))

  ;; resolve the polygon p
  (setf polygons (resolve p notch c-list))

  ;; for each polygon, recursively call acd and collect
  (loop for polygon in polygons
    do  (setf new-polygons (acd polygon tau))
        (loop for new-polygon in new-polygons
          do  (setf polygon-list (add new-polygon polygon-list))))

  (values polygon-list)))


;; run test on ACD algorithm 
(defun test ()
  ;; testing parameters
  (let* ( (output-file "output.txt")
          (polygons)
          (tau 1)
          (p (list  (list  1.0  1.0)
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
                    (list  2.0  1.0))))
  (setf polygons (acd p tau))
  (format t "~%>> Input  Number of Vertices: ~S" (length p))
  (format t "~%>> Output Number of Polygons: ~S" (length polygons))
  (format t "~%>> Writing Polygons to file:  ")
  (format t output-file)
  (with-open-file (str output-file
                      :direction :output
                      :if-exists :supersede
                      :if-does-not-exist :create)
    (loop for polygon in polygons
      do  (format str (write-to-string polygon))
          (format str "~%")))))

