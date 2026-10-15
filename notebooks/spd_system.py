import math

import pygame

pygame.init()

# ----------------------------
# Window
# ----------------------------
WIDTH, HEIGHT = 1000, 600
screen = pygame.display.set_mode((WIDTH, HEIGHT))
pygame.display.set_caption("Spring-Mass-Damper")

clock = pygame.time.Clock()
font = pygame.font.SysFont(None, 28)

# ----------------------------
# Physics
# ----------------------------
m = 1.0
k = 20.0
c = 2.0

force_amp = 0.0
force_freq = 1.0

x = 150.0
v = 0.0

equilibrium = 350

# ----------------------------
# Visual
# ----------------------------
wall_x = 80
mass_w = 80
mass_h = 60

dragging = False
paused = False

history = []

t = 0.0


# ----------------------------
# Helpers
# ----------------------------
def draw_text(text, x, y):
    img = font.render(text, True, (240, 240, 240))
    screen.blit(img, (x, y))


def draw_spring(x1, x2, y):
    coils = 14

    points = [(x1, y)]

    for i in range(1, coils):
        px = x1 + (x2 - x1) * i / coils

        py = y + (20 if i % 2 else -20)

        points.append((px, py))

    points.append((x2, y))

    pygame.draw.lines(screen, (100, 220, 255), False, points, 3)


# ----------------------------
# Main loop
# ----------------------------
running = True

while running:
    dt = clock.tick(60) / 1000
    t += dt

    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False

        elif event.type == pygame.KEYDOWN:
            if event.key == pygame.K_r:
                x = 150
                v = 0

            elif event.key == pygame.K_SPACE:
                paused = not paused

            elif event.key == pygame.K_UP:
                c += 0.5

            elif event.key == pygame.K_DOWN:
                c = max(0, c - 0.5)

            elif event.key == pygame.K_RIGHT:
                k += 2

            elif event.key == pygame.K_LEFT:
                k = max(1, k - 2)

            elif event.key == pygame.K_f:
                force_amp += 5

            elif event.key == pygame.K_g:
                force_amp = max(0, force_amp - 5)

        elif event.type == pygame.MOUSEBUTTONDOWN:
            mx, my = pygame.mouse.get_pos()

            mass_x = equilibrium + x

            rect = pygame.Rect(mass_x, HEIGHT // 2 - mass_h // 2, mass_w, mass_h)

            if rect.collidepoint(mx, my):
                dragging = True

        elif event.type == pygame.MOUSEBUTTONUP:
            dragging = False

    # ------------------------
    # Mouse drag
    # ------------------------
    if dragging:
        mx, _ = pygame.mouse.get_pos()

        x = mx - equilibrium
        v = 0

    # ------------------------
    # Physics update
    # ------------------------
    elif not paused:
        forcing = force_amp * math.sin(2 * math.pi * force_freq * t)

        a = (-k * x - c * v + forcing) / m

        v += a * dt
        x += v * dt

    # ------------------------
    # Graph history
    # ------------------------
    history.append(x)

    if len(history) > 300:
        history.pop(0)

    # ------------------------
    # Energies
    # ------------------------
    kinetic = 0.5 * m * v * v
    potential = 0.5 * k * x * x
    total_energy = kinetic + potential

    critical_c = 2 * math.sqrt(k * m)

    if c < critical_c:
        state = "UNDERDAMPED"
    elif abs(c - critical_c) < 0.1:
        state = "CRITICAL"
    else:
        state = "OVERDAMPED"

    # ------------------------
    # Draw
    # ------------------------
    screen.fill((30, 30, 35))

    # Wall
    pygame.draw.rect(screen, (170, 170, 170), (wall_x - 25, 80, 25, 300))

    mass_x = equilibrium + x
    mass_y = HEIGHT // 2 - mass_h // 2

    # Spring
    draw_spring(wall_x, mass_x, HEIGHT // 2)

    # Mass
    pygame.draw.rect(
        screen, (255, 180, 80), (mass_x, mass_y, mass_w, mass_h), border_radius=8
    )

    # Equilibrium marker
    pygame.draw.line(screen, (100, 255, 100), (equilibrium, 100), (equilibrium, 450), 2)

    # ------------------------
    # Graph
    # ------------------------
    graph_y = 500

    pygame.draw.line(screen, (120, 120, 120), (0, graph_y), (WIDTH, graph_y), 1)

    if len(history) > 1:
        pts = []

        for i, value in enumerate(history):
            gx = 600 + i

            gy = graph_y - value * 0.25

            pts.append((gx, gy))

        pygame.draw.lines(screen, (255, 100, 100), False, pts, 2)

    # ------------------------
    # Text
    # ------------------------
    draw_text(f"Mass M: {m:.1f}", 20, 20)
    draw_text(f"Spring K: {k:.1f}", 20, 50)
    draw_text(f"Damping C: {c:.1f}", 20, 80)

    draw_text(f"Displacement: {x:.2f}", 20, 130)
    draw_text(f"Velocity: {v:.2f}", 20, 160)

    draw_text(f"Energy: {total_energy:.2f}", 20, 190)

    draw_text(f"Force Amp: {force_amp:.1f}", 20, 220)

    draw_text(f"State: {state}", 20, 260)

    draw_text("Drag block with mouse", 20, 320)

    draw_text("UP/DOWN : damping", 20, 350)

    draw_text("LEFT/RIGHT : spring", 20, 380)

    draw_text("F/G : forcing", 20, 410)

    draw_text("SPACE : pause", 20, 440)

    draw_text("R : reset", 20, 470)

    pygame.display.flip()

pygame.quit()
