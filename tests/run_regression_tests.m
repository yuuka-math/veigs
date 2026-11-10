function report = run_regression_tests(repo_dir)
%RUN_REGRESSION_TESTS Check verified bounds against analytic eigenvalue oracles.
% Initialize INTLAB before calling. With no argument, test this repository.
% The returned report is a struct; no files are written by this test suite.
if nargin < 1
    repo_dir = fileparts(fileparts(mfilename('fullpath')));
end
original_path = path;
restore_path = onCleanup(@() path(original_path)); %#ok<NASGU>
addpath(repo_dir, '-begin');
clear veig veigs GetInertia
if isempty(which('intval'))
    error('veigs_tests:NoINTLAB', 'Initialize INTLAB before running this suite.');
end

cases = repmat(struct('name', '', 'passed', false, 'error', ''), 0, 1);
cases(end+1) = run_case('veig_scalar_exact', @scalar_case);
cases(end+1) = run_case('veig_dense_exact', @dense_case);
cases(end+1) = run_case('veig_generalized_exact', @generalized_case);
cases(end+1) = run_case('veig_repeated_eigenvalues', @repeated_case);
cases(end+1) = run_case('veig_interval_scalar_endpoints', @interval_scalar_case);
cases(end+1) = run_case('veig_negative_dense_extreme', @negative_dense_case);

selectors = {2, 'sa', 'la'};
for reflected = 0:1
    for k = 1:numel(selectors)
        name = sprintf('veigs_cluster_reflection_%d_selector_%d', reflected, k);
        cases(end+1) = run_case(name, @() cluster_case(reflected, selectors{k})); %#ok<AGROW>
    end
end
cases(end+1) = run_case('veigs_nonpositive_B_bound_rejected', @metric_guard_case);
cases(end+1) = run_case('veigs_distant_uncertainty_right', @() distant_uncertainty_case(1));
cases(end+1) = run_case('veigs_distant_uncertainty_left', @() distant_uncertainty_case(-1));
cases(end+1) = run_case('veigs_large_positive_extreme', @() large_extreme_case(1));
cases(end+1) = run_case('veigs_large_negative_extreme', @() large_extreme_case(-1));
cases(end+1) = run_case('veigs_nonfinite_rough_bound_rejected', @finite_guard_case);

report = struct('matlab_version', version, 'veigs_source', which('veigs'), ...
    'cases', cases, 'pass_count', sum([cases.passed]), ...
    'fail_count', sum(~[cases.passed]), 'all_passed', all([cases.passed]));
fprintf('%s\n', jsonencode(report));
if ~report.all_passed
    error('veigs_tests:RegressionFailure', '%d of %d regression checks failed.', ...
        report.fail_count, numel(cases));
end
end

function entry = run_case(name, check)
entry = struct('name', name, 'passed', false, 'error', '');
try
    check();
    entry.passed = true;
catch ME
    entry.error = sprintf('%s: %s', ME.identifier, ME.message);
end
end

function scalar_case()
[z, idx] = veig(7, 2);
assert_encloses(z, idx, 3.5, 3.5, 1);
end

function dense_case()
% Eigenvectors [1;-1] and [1;1] give eigenvalues 1 and 9 exactly.
[z, idx] = veig([5 4; 4 5], eye(2));
assert_encloses(z, idx, [1;9], [1;9], [1;2]);
end

function generalized_case()
% Shared eigenvectors give exact ratios 1/1 and 18/9.
[z, idx] = veig([9.5 8.5; 8.5 9.5], [5 4; 4 5]);
assert_encloses(z, idx, [1;2], [1;2], [1;2]);
end

function repeated_case()
[z, idx] = veig(diag([1 1 3 3]), eye(4), [2 3]);
assert_encloses(z, idx, [1;1;3;3], [1;1;3;3], [2;3]);
end

function interval_scalar_case()
% The exact eigenvalue range for this family is [7/4,9/2].
[z, idx] = veig(infsup(7, 9), infsup(2, 4));
assert_encloses(z, idx, 1.75, 4.5, 1);
end

function negative_dense_case()
% The smaller eigenvalue is -63-sqrt(4097). Certify endpoint positions
% using interval polynomial signs, without relying on approximate EIG.
[z, idx] = veig([-64 -64; -64 -62], eye(2), 1);
location = find(idx == 1, 1);
assert(~isempty(location), 'The requested eigenvalue index is missing.');
lo = inf(z(location)); hi = sup(z(location));
assert(isfinite(lo) && isfinite(hi), 'Nonfinite returned endpoint.');
qlo = (intval(lo)+63)^2-4097;
qhi = (intval(hi)+63)^2-4097;
assert(lo <= -63 && inf(qlo) >= 0, 'Lower endpoint misses the exact root.');
assert(hi >= -63 || sup(qhi) <= 0, 'Upper endpoint misses the exact root.');
end

function cluster_case(reflected, selector)
% The variable 2x2 block has eigenvalues (5 +/- sqrt(1+4*t^2))/2,
% |t| <= the represented double 0.1. These extrema enclose the whole family.
B = eye(3);
A = intval(diag([1,2,3]));
A(2,3) = infsup(-0.1,0.1); A(3,2) = A(2,3);
gap = sqrt(1+4*intval(0.1)^2);
oracle_lo = [1;inf((5-gap)/2);3];
oracle_hi = [1;2;sup((5+gap)/2)];
if reflected
    A = 4*intval(B)-A;
    reflected_oracle = 4-infsup(oracle_lo,oracle_hi);
    oracle_lo = flipud(inf(reflected_oracle));
    oracle_hi = flipud(sup(reflected_oracle));
end
[z, idx] = veigs(A, B, 2, selector);
assert(numel(idx) >= 2, 'Fewer than the requested two eigenvalues returned.');
assert_encloses(z, idx, oracle_lo, oracle_hi, []);
end

function metric_guard_case()
assert_unsafe_scale(@() veigs(diag([1,2,3]), 1e-300*eye(3), 1, 'sa'));
end

function distant_uncertainty_case(direction)
% Only a distant eigenvalue is uncertain. A dimension multiplier followed
% by binary halving skips a feasible residual radius between 1 and 1.2.
if direction > 0
    A = infsup(diag([0,10,11.2]), diag([2,10,11.2]));
    oracle_lo = [0;10;11.2]; oracle_hi = [2;10;11.2];
else
    A = infsup(diag([-11.2,-10,-2]), diag([-11.2,-10,0]));
    oracle_lo = [-11.2;-10;-2]; oracle_hi = [-11.2;-10;0];
end
[z, idx] = veigs(A, eye(3), 1, direction*9.9);
assert_encloses(z, idx, oracle_lo, oracle_hi, 2);
end

function large_extreme_case(direction)
% Entries and eigenvalues are finite, but 2*norm(A,inf) would overflow.
A = direction*pow2(1023)*eye(2);
% Both signs exercise the large finite global bound; the negative case
% also enters the lower rough search during local separation.
[z, idx] = veigs(A, 4*eye(2), 1, 'la');
exact = direction*pow2(1021)*ones(2,1);
assert_encloses(z, idx, exact, exact, 2);
end

function finite_guard_case()
% The infinity-norm quotient is unrepresentable. Reject it instead of
% allowing Inf <= Inf to keep a rough-bound search alive indefinitely.
assert_unsafe_scale(@() veigs(pow2(1023)*eye(2), 0.5*eye(2), 1, 'la'));
end

function assert_unsafe_scale(call)
try
    call();
catch ME
    assert(strcmp(ME.identifier, 'veigs:UnsafeScale'), ...
        'Expected veigs:UnsafeScale, received %s: %s', ME.identifier, ME.message);
    return
end
error('veigs_tests:MissingError', 'Expected veigs:UnsafeScale was not thrown.');
end

function assert_encloses(z, idx, oracle_lo, oracle_hi, requested)
lo = inf(z(:)); hi = sup(z(:)); idx = idx(:);
assert(numel(lo) == numel(idx), 'Bounds and index lengths differ.');
assert(all(isfinite(lo)) && all(isfinite(hi)), 'Nonfinite returned endpoint.');
assert(all(idx >= 1 & idx <= numel(oracle_lo) & idx == fix(idx)), ...
    'Invalid returned eigenvalue index.');
assert(numel(unique(idx)) == numel(idx), 'Duplicate eigenvalue indices.');
assert(all(ismember(requested, idx)), 'A requested eigenvalue index is missing.');
assert(all(lo <= oracle_lo(idx)) && all(hi >= oracle_hi(idx)), ...
    'Returned bounds miss an analytic eigenvalue endpoint.');
end
