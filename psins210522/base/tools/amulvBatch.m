function vo = amulvBatch(att, vi)
% vn transform by att to get vb i.e. vb = qmulv(qconj(a2qua(att)),vn). 
% See also a2qua, qconj, qmulv.
    att2 = att*0.5;
    s = sin(att2); c = cos(att2);
    sp = s(:,1); sr = s(:,2); sy = s(:,3); 
    cp = c(:,1); cr = c(:,2); cy = c(:,3); 
    q = [ cp.*cr.*cy - sp.*sr.*sy, sp.*cr.*cy - cp.*sr.*sy, cp.*sr.*cy + sp.*cr.*sy, cp.*cr.*sy + sp.*sr.*cy ];
    q = [q(:,1), -q(:,2), -q(:,3), -q(:,4)];
    qo1 =              - q(:,2) .* vi(:,1) - q(:,3) .* vi(:,2) - q(:,4) .* vi(:,3);
    qo2 = q(:,1) .* vi(:,1)                + q(:,3) .* vi(:,3) - q(:,4) .* vi(:,2);
    qo3 = q(:,1) .* vi(:,2)                + q(:,4) .* vi(:,1) - q(:,2) .* vi(:,3);
    qo4 = q(:,1) .* vi(:,3)                + q(:,2) .* vi(:,2) - q(:,3) .* vi(:,1);
    vo = vi;
    vo(:,1) = -qo1 .* q(:,2) + qo2 .* q(:,1) - qo3 .* q(:,4) + qo4 .* q(:,3);
    vo(:,2) = -qo1 .* q(:,3) + qo3 .* q(:,1) - qo4 .* q(:,2) + qo2 .* q(:,4);
    vo(:,3) = -qo1 .* q(:,4) + qo4 .* q(:,1) - qo2 .* q(:,3) + qo3 .* q(:,2);
