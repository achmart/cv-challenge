function [I1r, I2r, T1, T2, T3] = RectifyImages(I1, I2, K, R2, T2, R3, T3, do_plot)
    % This function rectifies two stereo images and returns a 
    % Translation for Derectifying a virtual viewpoint.
    
    % For computing the rectifying homographies it uses the a method
    % described by A.Fusiello et al. in 1999 
    % (http://www.diegm.uniud.it/fusiello/papers/00120016.pdf)
    
    % It takes to input images I1 & I2, the Kalibration matrix K, the 
    % Rotation R2 and Translation T2 between the camera frames and the
    % Rotation R3 and Translation T3 between the first camera frame and
    % the virtual viewpoint.
    % It returns the two rectified images I1r & I2r and the corresponding 
    % translation of the retinal plane T1 & T2 for the images and T3 for
    % the virtual viewpoint.
    
    % used parameters
    A1 = K;
    R1 = eye(3);
    t1 = zeros(3);
    Po1 = A1 * [R1, t1];
    %R2 = R2;
    t2 = T2;
    A2 = K;
    Po2 = A2 * [R2, t2];
    A3 = K;
    %R3 = R3;
    t3 = T3;
    Po3 = A3 * [R3,t3];

    % optical centers (unchanged)
    c1 = - R1'*(A1\Po1(:,4));
    c2 = - R2'*(A2\Po2(:,4));
    c3 = - R3'*(A3\Po3(:,4));

    % new x axis (baseline, from c1 to c2)
    v1 = (c2-c1);
    % new y axes (orthogonal to old z and new x)
    v2 = cross(R1(3,:)',v1);
    % new z axes (no choice, orthogonal to baseline and y)
    v3 = cross(v1,v2);

    % new extrinsic (translation unchanged)
    R = [v1'/norm(v1)
        v2'/norm(v2)
        v3'/norm(v3)];

    % new intrinsic (arbitrary) 
    An1 = A2;
    An1(1,2)=0;
    An2 = A2;
    An2(1,2)=0;
    An3 = A3;
    An3(1,2)=0;

    % translate image centers 
    An1(1,3)=An1(1,3);
    An1(2,3)=An1(2,3);
    An2(1,3)=An2(1,3);
    An2(2,3)=An2(2,3);
    An3(1,3)=An3(1,3);
    An3(2,3)=An3(2,3);

    % new projection matrices
    Pn1 = An1 * [R -R*c1 ];
    Pn2 = An2 * [R -R*c2 ];
    Pn3 = An3 * [R -R*c3 ];

    % rectifying image transformation
    T1 = Pn1(1:3,1:3)/(Po1(1:3,1:3));
    T2 = Pn2(1:3,1:3)/(Po2(1:3,1:3));
    T3 = Pn3(1:3,1:3)/(Po3(1:3,1:3));
    
    % find the smallest bounding box containining both images
    bb = mcbb(size(I1),size(I2), T1, T2);
    
    % warp RGB channels,
    parfor c = 1:size(I1, 3)

        % Warp LEFT
        [I1r(:,:,c),~,~] = imagewarp(I1(:,:,c), T1, 'bilinear', bb);

        % Warp RIGHT
        [I2r(:,:,c),~,~] = imagewarp(I2(:,:,c), T2, 'bilinear', bb);

    end
    
    %% Plot rectified images
    %if(do_plot)
    %    figure
    %    imshowpair(I1r,I2r);
    %end
end