//! The section of a face on a surface of revolution by a plane through its axis, in closed
//! form: the kernel's [`solidcore::geom::axial_section::axial_section`] on records' numbers
//! (#17496, #17498, #17499, #17500). Ops converts.

use solidcore::geom::axial_section::{
    self as k, AxialPiece as KPiece, CommonAxis as KCommonAxis, Revolved as KRevolved,
};
use solidcore::geom::vector3d::Vec3;

/// A surface of revolution ([`k::Revolved`]).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Revolved {
    /// A cylinder: a point on its axis, the axis and the radius.
    Cylinder {
        /// A point on the axis.
        origin: [f64; 3],
        /// The axis direction (any length).
        axis: [f64; 3],
        /// The radius.
        radius: f64,
    },
    /// A cone: its apex, its axis pointing into the cone and its half angle (radians).
    Cone {
        /// The apex.
        apex: [f64; 3],
        /// The axis direction from the apex into the cone (any length).
        axis: [f64; 3],
        /// The angle between the axis and each line of the cone (radians).
        half_angle: f64,
    },
    /// A sphere.
    Sphere {
        /// The centre.
        centre: [f64; 3],
        /// The radius.
        radius: f64,
    },
    /// A plane: revolved about any axis along its normal.
    Plane {
        /// A point on the plane.
        point: [f64; 3],
        /// Its normal (any length).
        normal: [f64; 3],
    },
    /// A torus: its centre, axis, ring (major) and tube (minor) radii.
    Torus {
        /// The centre of the ring.
        centre: [f64; 3],
        /// The axis direction (any length).
        axis: [f64; 3],
        /// The radius of the circle the tube's centre runs on.
        major: f64,
        /// The tube's radius.
        minor: f64,
    },
}

/// A piece of an axial section, in the plane ([`k::AxialPiece`]).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum AxialPiece {
    /// A segment of a line of a cylinder or a cone.
    Line {
        /// Where it starts.
        start: [f64; 3],
        /// Where it ends.
        end: [f64; 3],
    },
    /// An arc, counter-clockwise about the plane's normal from `start` to `end`.
    Arc {
        /// The circle's centre.
        centre: [f64; 3],
        /// Where it starts.
        start: [f64; 3],
        /// Where it ends.
        end: [f64; 3],
    },
    /// A whole circle.
    Circle {
        /// The circle's centre.
        centre: [f64; 3],
        /// Its radius.
        radius: f64,
    },
}

fn vec([x, y, z]: [f64; 3]) -> Vec3 {
    Vec3::new(x, y, z)
}

fn surface(s: &Revolved) -> KRevolved {
    match *s {
        Revolved::Cylinder {
            origin,
            axis,
            radius,
        } => KRevolved::Cylinder {
            origin: vec(origin),
            axis: vec(axis),
            radius,
        },
        Revolved::Cone {
            apex,
            axis,
            half_angle,
        } => KRevolved::Cone {
            apex: vec(apex),
            axis: vec(axis),
            half_angle,
        },
        Revolved::Sphere { centre, radius } => KRevolved::Sphere {
            centre: vec(centre),
            radius,
        },
        Revolved::Plane { point, normal } => KRevolved::Plane {
            point: vec(point),
            normal: vec(normal),
        },
        Revolved::Torus {
            centre,
            axis,
            major,
            minor,
        } => KRevolved::Torus {
            centre: vec(centre),
            axis: vec(axis),
            major,
            minor,
        },
    }
}

fn piece(p: KPiece) -> AxialPiece {
    match p {
        KPiece::Line { start, end } => AxialPiece::Line {
            start: start.to_array(),
            end: end.to_array(),
        },
        KPiece::Arc { centre, start, end } => AxialPiece::Arc {
            centre: centre.to_array(),
            start: start.to_array(),
            end: end.to_array(),
        },
        KPiece::Circle { centre, radius } => AxialPiece::Circle {
            centre: centre.to_array(),
            radius,
        },
    }
}

/// How far a line leaves the plane through `plane_origin` spanned by `plane_x` and `plane_y`:
/// `(tilt, offset)`, the line's unit direction's component along the plane's unit normal, and the height of `line_point` above the plane (the kernel's
/// [`k::line_plane_deviation`]). `None` for a number that is not finite, axes that do not span a
/// plane or a direction of no length.
#[must_use]
pub fn line_plane_deviation(
    plane_origin: [f64; 3],
    plane_x: [f64; 3],
    plane_y: [f64; 3],
    line_point: [f64; 3],
    line_dir: [f64; 3],
) -> Option<(f64, f64)> {
    k::line_plane_deviation(
        vec(plane_origin),
        vec(plane_x),
        vec(plane_y),
        vec(line_point),
        vec(line_dir),
    )
}

/// The range `(lowest, highest)` of the speed across the plane spanned by `plane_x` and
/// `plane_y` of `points` on it, under the screw motion about the line through `axis_point` along
/// `axis_dir` rising `pitch` per turn (the kernel's [`k::screw_speed_range`]): a section on the
/// plane sweeps without passing back over itself while the speed keeps one sign over it. `None`
/// for a number that is not finite, axes that do not span a plane, a direction of no length or
/// no points.
#[must_use]
pub fn screw_speed_range(
    plane_x: [f64; 3],
    plane_y: [f64; 3],
    axis_point: [f64; 3],
    axis_dir: [f64; 3],
    pitch: f64,
    points: &[[f64; 3]],
) -> Option<(f64, f64)> {
    let points: Vec<Vec3> = points.iter().map(|&p| vec(p)).collect();
    k::screw_speed_range(
        vec(plane_x),
        vec(plane_y),
        vec(axis_point),
        vec(axis_dir),
        pitch,
        &points,
    )
}

/// Whether the plane through `plane_point` square to `plane_normal` holds the axis of `on`
/// within `tol` (the model's tolerance) wherever it passes `extent` (a box around the face), or
/// a sphere's centre (the kernel's [`k::holds_axis`]): whether [`axial_section`] sections the
/// face there.
#[must_use]
pub fn holds_axis(
    on: &Revolved,
    plane_point: [f64; 3],
    plane_normal: [f64; 3],
    extent: &[[f64; 3]],
    tol: f64,
) -> bool {
    let extent: Vec<Vec3> = extent.iter().map(|p| vec(*p)).collect();
    k::holds_axis(
        &surface(on),
        vec(plane_point),
        vec(plane_normal),
        &extent,
        tol,
    )
}

/// The section of a face on `surface` by the plane through `plane_point` square to
/// `plane_normal`, when the plane holds the surface's axis within `tol` (the model's
/// tolerance) wherever it passes `extent` (a box around the face, [`holds_axis`]); `cuts` are where the face's boundary crosses the
/// plane, and `on_face` says whether a point of the surface is the face's (the kernel's
/// [`k::axial_section`]). `Ok(None)` when the plane does not hold the axis or a number is not
/// usable; `Ok(Some(vec![]))` when the plane misses the face.
///
/// # Errors
///
/// `on_face`'s error.
pub fn axial_section<E>(
    on: &Revolved,
    plane_point: [f64; 3],
    plane_normal: [f64; 3],
    extent: &[[f64; 3]],
    tol: f64,
    cuts: &[[f64; 3]],
    mut on_face: impl FnMut([f64; 3]) -> Result<bool, E>,
) -> Result<Option<Vec<AxialPiece>>, E> {
    let cuts: Vec<Vec3> = cuts.iter().map(|p| vec(*p)).collect();
    let extent: Vec<Vec3> = extent.iter().map(|p| vec(*p)).collect();
    let made = k::axial_section(
        &surface(on),
        vec(plane_point),
        vec(plane_normal),
        &extent,
        tol,
        &cuts,
        |p| on_face(p.to_array()),
    )?;
    Ok(made.map(|pieces| pieces.into_iter().map(piece).collect()))
}

/// The axis a set of faces share ([`k::CommonAxis`] in records' numbers), as [`common_axis`]
/// finds it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum CommonAxis {
    /// Every face is revolved about the line through `origin` along the unit `axis`; `by` is the
    /// face whose axis it is (an index into the faces passed).
    About {
        /// A point on the axis.
        origin: [f64; 3],
        /// The unit axis.
        axis: [f64; 3],
        /// The face that gave it.
        by: usize,
    },
    /// No face gives an axis: only planes and spheres.
    Undetermined,
    /// Face `face` is not revolved about the first axis tried, and no other face's axis holds
    /// every face.
    NotAbout {
        /// The face, as an index into the faces passed.
        face: usize,
    },
    /// Face `face` is no surface of revolution or plane (`None` passed for it).
    Unproved {
        /// The face, as an index into the faces passed.
        face: usize,
    },
}

/// The axis every face is revolved about within `tol` (the kernel's [`k::common_axis`]), each
/// face given by its surface (`None` for one that is no surface of revolution or plane) and the
/// corners of a box around it.
#[must_use]
pub fn common_axis(faces: &[(Option<Revolved>, Vec<[f64; 3]>)], tol: f64) -> CommonAxis {
    let mut surfaces = Vec::with_capacity(faces.len());
    for (face, (on, extent)) in faces.iter().enumerate() {
        let Some(on) = on else {
            return CommonAxis::Unproved { face };
        };
        surfaces.push((surface(on), extent.iter().map(|p| vec(*p)).collect()));
    }
    match k::common_axis(&surfaces, tol) {
        KCommonAxis::About { origin, axis, by } => CommonAxis::About {
            origin: origin.to_array(),
            axis: axis.to_array(),
            by,
        },
        KCommonAxis::Undetermined => CommonAxis::Undetermined,
        KCommonAxis::NotAbout { surface } => CommonAxis::NotAbout { face: surface },
    }
}
