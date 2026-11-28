//! The axis a set of surfaces of revolution share (`geom::axial_section::common_axis`), and a
//! plane's section by a plane its normal lies in (`Revolved::Plane`): two tubes' cylinders,
//! bevel cone, land plane, a torus and a sphere on the axis share it; a cylinder moved off it, or
//! a plane leaning across it, is named; planes and spheres alone give none; a land's section by a
//! plane through the axis is its two lines, kept where the face is (T49, #18500).
#![allow(clippy::unwrap_used, clippy::panic, clippy::indexing_slicing)] // test helpers
use solidcore::geom::axial_section::{
    axial_section, common_axis, AxialPiece, CommonAxis, Revolved,
};
use solidcore::geom::vector3d::Vec3;

const TOL: f64 = 1e-5;

fn v(x: f64, y: f64, z: f64) -> Vec3 {
    Vec3::new(x, y, z)
}

/// The corners of the box from (-30, -48, -48) to (30, 48, 48): two tubes end to end.
fn tube_box() -> Vec<Vec3> {
    let mut out = Vec::new();
    for x in [-30.0, 30.0] {
        for y in [-48.0, 48.0] {
            for z in [-48.0, 48.0] {
                out.push(v(x, y, z));
            }
        }
    }
    out
}

fn cylinder(y: f64, radius: f64) -> Revolved {
    Revolved::Cylinder {
        origin: v(-30.0, y, 0.0),
        axis: v(1.0, 0.0, 0.0),
        radius,
    }
}

fn land(normal: Vec3) -> Revolved {
    Revolved::Plane {
        point: v(-1.0, 0.0, 0.0),
        normal,
    }
}

fn shared() -> Vec<Revolved> {
    vec![
        cylinder(0.0, 48.0),
        Revolved::Cone {
            apex: v(23.248_711_305_964_28, 0.0, 0.0),
            axis: v(-1.0, 0.0, 0.0),
            half_angle: 60.0_f64.to_radians(),
        },
        land(v(1.0, 0.0, 0.0)),
        cylinder(0.0, 40.0),
        Revolved::Torus {
            centre: v(5.0, 0.0, 0.0),
            axis: v(-1.0, 0.0, 0.0),
            major: 44.0,
            minor: 3.0,
        },
        Revolved::Sphere {
            centre: v(12.0, 0.0, 0.0),
            radius: 10.0,
        },
    ]
}

fn boxed(surfaces: Vec<Revolved>) -> Vec<(Revolved, Vec<Vec3>)> {
    surfaces.into_iter().map(|s| (s, tube_box())).collect()
}

#[test]
fn surfaces_about_one_axis_share_it() {
    match common_axis(&boxed(shared()), TOL) {
        CommonAxis::About { origin, axis, by } => {
            assert_eq!(by, 0);
            assert!(
                origin.y.abs() < 1e-12 && origin.z.abs() < 1e-12,
                "{origin:?}"
            );
            assert!((axis.x.abs() - 1.0).abs() < 1e-12, "{axis:?}");
        }
        other => panic!("{other:?}"),
    }
    // the order of the surfaces changes which gives the axis, not the answer
    let mut turned = shared();
    turned.rotate_left(2);
    assert!(matches!(
        common_axis(&boxed(turned), TOL),
        CommonAxis::About { by: 1, .. }
    ));
}

#[test]
fn a_surface_off_the_axis_is_named() {
    // the bore moved 0.01 off the axis
    let mut off = shared();
    off[3] = cylinder(0.01, 40.0);
    assert_eq!(
        common_axis(&boxed(off), TOL),
        CommonAxis::NotAbout { surface: 3 }
    );
    // within the tolerance it is on it
    let mut near = shared();
    near[3] = cylinder(0.5 * TOL, 40.0);
    assert!(matches!(
        common_axis(&boxed(near), TOL),
        CommonAxis::About { .. }
    ));
    // the land leaning 0.001 rad across the axis: 0.06 off over the box's reach
    let mut leaning = shared();
    leaning[2] = land(v(1.0, 0.001, 0.0));
    assert_eq!(
        common_axis(&boxed(leaning), TOL),
        CommonAxis::NotAbout { surface: 2 }
    );
    // a sphere whose centre is off the axis
    let mut ball = shared();
    ball[5] = Revolved::Sphere {
        centre: v(12.0, 0.0, 0.01),
        radius: 10.0,
    };
    assert_eq!(
        common_axis(&boxed(ball), TOL),
        CommonAxis::NotAbout { surface: 5 }
    );
}

#[test]
fn planes_and_spheres_alone_give_no_axis() {
    let only = vec![
        land(v(1.0, 0.0, 0.0)),
        Revolved::Sphere {
            centre: v(12.0, 0.0, 0.0),
            radius: 10.0,
        },
    ];
    assert_eq!(common_axis(&boxed(only), TOL), CommonAxis::Undetermined);
    assert_eq!(common_axis(&[], TOL), CommonAxis::Undetermined);
}

#[test]
fn a_land_cut_through_the_axis_is_its_two_lines() {
    // the land x = -1 between the bore (r 40) and the bevel's foot (r 42), cut by z = 0
    let cuts = [
        v(-1.0, 40.0, 0.0),
        v(-1.0, 42.0, 0.0),
        v(-1.0, -40.0, 0.0),
        v(-1.0, -42.0, 0.0),
    ];
    let on_face = |p: Vec3| -> Result<bool, ()> {
        let r = p.y.hypot(p.z);
        Ok((40.0..=42.0).contains(&r) && (p.x + 1.0).abs() < 1e-9)
    };
    let pieces = axial_section(
        &land(v(1.0, 0.0, 0.0)),
        v(0.0, 0.0, 0.0),
        v(0.0, 0.0, 1.0),
        &tube_box(),
        TOL,
        &cuts,
        on_face,
    )
    .unwrap()
    .unwrap();
    let mut ys: Vec<[f64; 2]> = pieces
        .iter()
        .map(|p| match *p {
            AxialPiece::Line { start, end } => {
                assert!(
                    (start.x + 1.0).abs() < 1e-12 && (end.x + 1.0).abs() < 1e-12,
                    "{p:?}"
                );
                [start.y.min(end.y), start.y.max(end.y)]
            }
            _ => panic!("{p:?}"),
        })
        .collect();
    ys.sort_by(|a, b| a[0].total_cmp(&b[0]));
    assert_eq!(ys, vec![[-42.0, -40.0], [40.0, 42.0]]);
    // a cutting plane the land's normal leaves gives none
    let none = axial_section(
        &land(v(1.0, 0.0, 0.0)),
        v(0.0, 0.0, 0.0),
        v(1.0, 0.0, 0.0),
        &tube_box(),
        TOL,
        &cuts,
        on_face,
    );
    assert_eq!(none.unwrap(), None);
}
