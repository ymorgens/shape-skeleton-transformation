function neglogaffinity = affinity(a,a_normal,p,p_normal,desired_vonMises_b,logconstant)
% neglogaffinity = affinity(a,a_normal,p,p_normal,desired_vonMises_b,logconstant)
% 
% Angular affinity between a(axis point/skeleton point) an p(shape point/boundary contour/shape polygon point); 
% shape normals assumed to point outward, so shapes should be drawn CLOCKWISE.

global vonMises_b;
if nargin<5
    desired_vonMises_b = vonMises_b;   
end

if nargin<6
    logconstant = log(2*pi*besseli(0,desired_vonMises_b));
end
    
p_to_a_penalty = 10000;
a_to_p_penalty = 10000;

a_to_p = p - a;

length_array = size(a,1);
if length_array > 1
    %---------------------------------------
    % cos_angular_deviation at shape point
    %---------------------------------------
    cos_angular_deviation_at_p = zeros(length_array,1);
    index1 = v_norm(p_normal) > 0;
    cos_angular_deviation_at_p = ( v_dot(-a_to_p,-p_normal) )./( v_norm(-a_to_p).*v_norm(-p_normal) );

    index2 = cos_angular_deviation_at_p < 0;
    cos_angular_deviation_at_p(index2) = p_to_a_penalty * cos_angular_deviation_at_p(index2);
    
    %---------------------------------------
    % cos_angular_deviation at skeletal point
    %---------------------------------------
    cos_angular_deviation_at_a = zeros(length_array,1);
    index3 = v_norm(a_normal) > 0;
    cos_angular_deviation_at_a = ( v_dot(a_to_p,a_normal) )./( v_norm(a_to_p).*v_norm(a_normal) );
    
    index4 = cos_angular_deviation_at_a < 0;
    cos_angular_deviation_at_a(index4) = a_to_p_penalty * cos_angular_deviation_at_a(index4);
    
    p_weighting = ones(length_array,1);
    % logconstant = log(2*pi*besseli(0,desired_vonMises_b)); %!!!! 2011-07-04, It takes too long time so I made it as parameter..from the upper level
    neglogaffinity = -desired_vonMises_b * (p_weighting .* cos_angular_deviation_at_p + cos_angular_deviation_at_a) ...
        + 2 * logconstant*ones(length_array,1);

else

    if norm(p_normal) > 0
        cos_angular_deviation_at_p = ...
            ( (-a_to_p(1)).*(-p_normal(1)) + (-a_to_p(2)).*(-p_normal(2)) ) / (norm(-a_to_p)*norm(-p_normal));
        %     dot(-a_to_p,-p_normal)/(norm(-a_to_p)*norm(-p_normal));

        % we use -p_normal to actually use the one that points INward, and
        % minus a_to_p (ie, "p to a"), so we comparing the two vectors that emerge from p, in towards a

        % Kluge:

        if cos_angular_deviation_at_p < 0 % angle is greater than 90deg
            %180 * acos(cos_angular_deviation_at_p)/pi
            cos_angular_deviation_at_p = p_to_a_penalty * cos_angular_deviation_at_p;
        end

    else
        cos_angular_deviation_at_p = 0;
    end

%---------------------------------------
% cos_angular_deviation at skeletal point
%---------------------------------------
if norm(a_normal) > 0
    cos_angular_deviation_at_a = ...
        ( (a_to_p(1)).*(a_normal(1)) + (a_to_p(2)).*(a_normal(2)) )/(norm(a_to_p)*norm(a_normal));
%     dot(a_to_p,a_normal)/(norm(a_to_p)*norm(a_normal));
    
    if cos_angular_deviation_at_a < 0
        cos_angular_deviation_at_a = a_to_p_penalty * cos_angular_deviation_at_a;
    end

else
    cos_angular_deviation_at_a = 0;
%     reverse_cos_angular_deviation_at_a = 0;
end


p_weighting = 1;
% logconstant = log(2*pi*besseli(0,desired_vonMises_b)); %!!!! 2011-07-04, It takes too long time so I made it as parameter..from the upper level
 neglogaffinity = -desired_vonMises_b * (p_weighting * cos_angular_deviation_at_p + cos_angular_deviation_at_a) + 2 * logconstant;

end

end


%=========
% v_norm
%=========
function r = v_norm(a)
% norm of [n*2] array
r = sqrt(a(:,1).^2 + a(:,2).^2);
end

%========
% v_dot
%========
function r = v_dot(a,b)
% dot of [n*2] arrays
r = (a(:,1)).*(b(:,1)) + (a(:,2)).*(b(:,2));
end

