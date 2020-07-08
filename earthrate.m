function omega_ie_n = earthrate(lat)
% earthrate: turn rate of the Earth in the navigation frame.
  
    omega_ie_n = (7.2921155e-5) .* [ (cos(lat));  0 ; -sin(lat) ]; 
end
