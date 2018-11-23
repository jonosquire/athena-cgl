function V = readHDF5(hdf5file)

F=hdf5info(hdf5file);

V.coord_sys = F.GroupHierarchy.Attributes(3).Value.Data;
V.t  = F.GroupHierarchy.Attributes(2).Value;
V.num_cycle = F.GroupHierarchy.Attributes(1).Value;
V.num_meshblocks = F.GroupHierarchy.Attributes(8).Value;
V.max_level = F.GroupHierarchy.Attributes(10).Value;
V.grid_num = F.GroupHierarchy.Attributes(7).Value.';
V.mesh_block_size = F.GroupHierarchy.Attributes(9).Value.';
names_struct = F.GroupHierarchy.Attributes(13).Value;
V.num_fields = length(names_struct);
dset_names = F.GroupHierarchy.Attributes(12).Value;
V.num_variables = F.GroupHierarchy.Attributes(11).Value.';
V.domain_bounds = [F.GroupHierarchy.Attributes(4).Value(1:2).';...
    F.GroupHierarchy.Attributes(5).Value(1:2).';...
    F.GroupHierarchy.Attributes(6).Value(1:2).'];


if V.max_level > 0
    error('Script not set up for mesh-refined simulations')
end
if ~strcmp(V.coord_sys,'cartesian')
    warning(['May not work with non-Cartesian coordinates'])
end

% Make some matrices to put things in
for fff = 1:V.num_fields
    names{fff} = names_struct(fff).data;
    V.(names{fff}) = zeros(V.grid_num);
end


% Put the data in the root grids
blocks = hdf5read(hdf5file,'/LogicalLocations');
mbs = int64(V.mesh_block_size);
varnum = 1;current_dset=1;
for fff = 1:V.num_fields
    if varnum == 1  % If you've finished the variable loop, reread the data
        current_data = hdf5read(hdf5file,['/' dset_names(current_dset).data]);
    end
    for mmm=1:V.num_meshblocks % Meshblocks get put at the 'LogicalLocations' 
        off = blocks(:,mmm).';
        ind_s = mbs.*off+1;
        ind_e = mbs.*off+mbs;
        V.(names{fff})(ind_s(1):ind_e(1), ind_s(2):ind_e(2), ind_s(3):ind_e(3)) = ...
            current_data(:,:,:,mmm,varnum);
    end
    varnum = varnum+1;
    if varnum>V.num_variables(current_dset)
        current_dset = current_dset + 1;
        varnum = 1;
    end
    
    
end


% Construct x, y and z vectors from individual meshblocks
% x vector 
blks = blocks(1,:);
current_blk = -1;
V.x = [];
coord_full = hdf5read(hdf5file,'/x1v');
for bbb = 1:length(blks)
    if blks(bbb)==current_blk+1
        V.x = [V.x; coord_full(:,bbb)];
        current_blk = current_blk + 1;
    end
end
% y vector 
blks = blocks(2,:);
current_blk = -1;
V.y = [];
coord_full = hdf5read(hdf5file,'/x2v');
for bbb = 1:length(blks)
    if blks(bbb)==current_blk+1
        V.y = [V.y; coord_full(:,bbb)];
        current_blk = current_blk + 1;
    end
end
% z vector 
blks = blocks(3,:);
current_blk = -1;
V.z = [];
coord_full = hdf5read(hdf5file,'/x3v');
for bbb = 1:length(blks)
    if blks(bbb)==current_blk+1
        V.z = [V.z; coord_full(:,bbb)];
        current_blk = current_blk + 1;
    end
end



end